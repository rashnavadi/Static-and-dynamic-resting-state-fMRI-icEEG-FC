%% written by Tahereh Rashnavadi
% Jan 2026


%% 04_compute_fc_fmri_ICE033.m
% Computes static and dynamic fMRI FC from bipolar seed time series after
% removing the shared monopolar electrode.
% Dynamic FC uses sliding windows and applies a forward time shift (lag)
% to account for hemodynamic delay (implemented as trimming the first lagTR samples).

% clear; clc;
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

% ---------- Load config + segTable ----------
if ~exist("C","var") || ~isfield(C,"subjectID") || strlength(C.subjectID)==0
    error("C.subjectID must be set by the calling script");
end

cfgPath = fullfile(C.projectRoot, "FC_rebuild", C.subjectID, sprintf("config_%s.mat", C.subjectID));
load(cfgPath, "C");

segPath = fullfile(C.outDir, sprintf("segTable_%s.mat", C.subjectID));
load(segPath, "segTable");


% ---------- Load EEG segments ----------
eegSegPath = fullfile(C.outDir, sprintf("EEG_segments_%s.mat", C.subjectID));
load(eegSegPath, "EEG");


% ---------- Load fMRI seeds ----------
fmriPath = fullfile(C.outDir, sprintf("FMRI_seedTS_%s.mat", C.subjectID));
load(fmriPath, "FMRI");

% ---------- Helpers ----------
addpath(C.codeDir);
addpath(fullfile(C.codeDir, "helpers"));

OUT = struct();

winTR  = round(C.winSec / C.TR);                     % 120/1.5=80
stepTR = round((C.winSec - C.overlapSec)/C.TR);      % 60/1.5=40
lagTR  = round(C.fmriLagSec / C.TR);                 % e.g., 6s/1.5=4

for r = 1:height(segTable)
    runID  = string(segTable.RunID(r));
    runKey = matlab.lang.makeValidName(runID);

    X = FMRI.(runKey).seedTS;           % [T x N]
    names = FMRI.(runKey).seedNames;
    names = upper(strrep(names, "_", "-"));
    names = canonical_bipolar_name(names);

    % Build bipolar names from EEG monopolar labels
    eeg_mono = EEG.(runKey).data_seg;
    monoLabels = string(EEG.(runKey).labels);

    [~, eegBipNames] = build_bipolar_signals_from_monopolar(eeg_mono, monoLabels);

    eegBipNames = upper(strrep(eegBipNames, "_", "-"));
    eegBipNames = canonical_bipolar_name(eegBipNames);

    % Keep only seeds that exist in EEG bipolars:
    keep = ismember(names, eegBipNames);

    N_before = numel(names);
    N_keep   = nnz(keep);
    N_drop   = N_before - N_keep;
    pct_keep = 100 * (N_keep / max(N_before,1));

    missingSeeds = names(~keep);

    fprintf("[MATCH] %s %s fMRI seeds: %d total | %d matched (%.1f%%) | %d dropped\n", ...
        C.subjectID, runID, N_before, N_keep, pct_keep, N_drop);

    if N_drop > 0
        fprintf("[INFO] %s %s dropping %d fMRI seeds not present in EEG:\n", ...
            C.subjectID, runID, N_drop);
        for m = 1:numel(missingSeeds)
            fprintf("       - %s\n", missingSeeds(m));
        end
    end

    % Now actually filter
    names = names(keep);
    X     = X(:, keep);
    N     = numel(names);

    % ---------- Seed-level de-overlap: keep only bipolars with unique contacts ----------
    % Example: keep DRFOP2-DRFOP3 and DRFOP4-DRFOP5, drop DRFOP3-DRFOP4 (shared contact)
    [seedKeep, droppedSeeds] = keep_unique_contact_bipolars(names);

    fprintf("[SEED-FILTER] %s %s kept %d/%d seeds (dropped %d due to shared contacts)\n", ...
        C.subjectID, runID, nnz(seedKeep), numel(names), numel(droppedSeeds));

    if ~isempty(droppedSeeds)
        for k = 1:numel(droppedSeeds)
            fprintf("       - %s\n", droppedSeeds(k));
        end
    end

    names = names(seedKeep);
    X     = X(:, seedKeep);
    N     = numel(names);

    if N < 2
        warning("[WARN] %s %s has <2 seeds after unique-contact filtering. Skipping run.", C.subjectID, runID);
        continue;
    end


    % ---------- Exclude shared-contact bipolar pairs ----------
    % Build a mask over the upper-triangle edges: keep only edges whose two bipolars
    % do NOT share any electrode contact (e.g., A-B with B-C is excluded).

    pairKeepMat = true(N,N);   % pairKeepMat(i,j)=true means keep edge i-j

    % Parse contacts for each bipolar name "A-B"
    c1 = strings(N,1);
    c2 = strings(N,1);
    for i = 1:N
        parts = split(names(i), "-");
        c1(i) = parts(1);
        c2(i) = parts(2);
    end

    % Mark shared-contact edges as false
    for i = 1:N
        for j = i+1:N
            shares = (c1(i)==c1(j)) || (c1(i)==c2(j)) || (c2(i)==c1(j)) || (c2(i)==c2(j));
            if shares
                pairKeepMat(i,j) = false;
                pairKeepMat(j,i) = false;
            end
        end
    end

    % Convert to vector mask matching uppertri_vec() order
    pairKeepVec = uppertri_vec(pairKeepMat);   % logical [nPairs x 1]
    nPairsKeep = nnz(pairKeepVec);
    fprintf("[FILTER] %s %s removing %d shared-contact edges (keeping %d)\n", ...
        C.subjectID, runID, (numel(pairKeepVec)-nPairsKeep), nPairsKeep);

    % ---------- Static FC ----------
    FC_static = corr(X, "Rows","complete");         % [N x N]
    FC_static = max(min(FC_static, 0.999999), -0.999999);
    FC_static_z = atanh(FC_static);                 % Fisher z

    % Also store filtered static edge vectors for downstream export/LMM
    staticVec   = uppertri_vec(FC_static);     % r
    staticVec   = staticVec(pairKeepVec);

    staticVec_z = uppertri_vec(FC_static_z);   % z
    staticVec_z = staticVec_z(pairKeepVec);

    % ---------- Dynamic FC ----------
    % Apply lag by shifting time windows forward on fMRI side:
    % We implement as analyzing X(lagTR+1:end,:) and computing windows there.
    if lagTR > 0
        Xlag = X(lagTR+1:end, :);
    else
        Xlag = X;
    end

    Tlag = size(Xlag,1);

    if Tlag < winTR
        warning("Run %s too short after lag (Tlag=%d < winTR=%d). Skipping dynamic FC.", runID, Tlag, winTR);
        winStarts = [];
        winEnds   = [];
        dynVec    = [];
        dynVec_z  = [];
        nWin      = 0;
    else
        [winStarts, winEnds] = sliding_window_indices(Tlag, winTR, stepTR, 0);
        nWin = numel(winStarts);
        dynVec   = nan(nPairsKeep, nWin);
        dynVec_z = nan(nPairsKeep, nWin);

        for w = 1:nWin
            xs = Xlag(winStarts(w):winEnds(w), :);
            Cw = corr(xs, "Rows","complete");
            v  = uppertri_vec(Cw);
            v  = v(pairKeepVec);   % <-- apply shared-contact exclusion here

            % safe Fisher-z
            v = max(min(v, 0.999999), -0.999999);

            dynVec(:,w)   = v;
            dynVec_z(:,w) = atanh(v);   % Fisher z per window (useful later)
        end
    end

    OUT.(runKey).seedNames    = names;
    OUT.(runKey).FC_static    = FC_static;
    OUT.(runKey).FC_static_z  = FC_static_z;
    OUT.(runKey).dynVec       = dynVec;
    OUT.(runKey).dynVec_z     = dynVec_z;
    OUT.(runKey).winStarts    = winStarts;
    OUT.(runKey).winEnds      = winEnds;
    OUT.(runKey).winStarts_orig = winStarts + lagTR;
    OUT.(runKey).winEnds_orig   = winEnds   + lagTR;
    OUT.(runKey).winTR        = winTR;
    OUT.(runKey).stepTR       = stepTR;
    OUT.(runKey).lagTR        = lagTR;
    % after remvoing the shared contacts in bipolar settings, e.g., A-B <-> B-C
    OUT.(runKey).pairKeepVec   = pairKeepVec;
    OUT.(runKey).staticVec     = staticVec;
    OUT.(runKey).staticVec_z   = staticVec_z;
    OUT.(runKey).nPairsKeep    = nPairsKeep;


    fprintf("[OK] fMRI FC %s %s: N=%d, windows=%d\n", C.subjectID, runID, N, nWin);
end

save(fullfile(C.outDir, sprintf("FMRI_FC_%s.mat", C.subjectID)), "OUT", "-v7.3");
fprintf("[OK] Saved %s\n", fullfile(C.outDir, sprintf("FMRI_FC_%s.mat", C.subjectID)));

function [keepMask, dropped] = keep_unique_contact_bipolars(names)
% keep_unique_contact_bipolars
% Keeps a subset of bipolar channels such that each monopolar contact
% is used at most once across the retained bipolars.
%
% Greedy strategy in the given order:
%   - keep a bipolar A-B if neither A nor B has been used already
%   - otherwise drop it
%
% Inputs:
%   names : string array, bipolar names like "DRFOP2-DRFOP3"
%
% Outputs:
%   keepMask : logical mask same length as names
%   dropped  : string array of dropped bipolar names

names = string(names(:));
n = numel(names);

keepMask = false(n,1);
dropped  = strings(0,1);

used = containers.Map('KeyType','char','ValueType','logical');

for i = 1:n
    nm = names(i);

    parts = split(nm, "-");
    if numel(parts) < 2
        % If formatting is unexpected, keep it (or change to drop if you prefer)
        keepMask(i) = true;
        continue;
    end

    a = char(strtrim(parts(1)));
    b = char(strtrim(parts(2)));

    aUsed = isKey(used, a);
    bUsed = isKey(used, b);

    if ~(aUsed || bUsed)
        keepMask(i) = true;
        used(a) = true;
        used(b) = true;
    else
        dropped(end+1,1) = nm; %#ok<AGROW>
    end
end
end

