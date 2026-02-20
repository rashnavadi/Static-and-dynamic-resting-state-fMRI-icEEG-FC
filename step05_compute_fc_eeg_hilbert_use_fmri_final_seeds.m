%% step05_compute_fc_eeg_hilbert_USE_FMRI_FINAL_SEEDS.m
% Written by Tahereh Rashnavadi
% Feb 2026
%
% Computes static + dynamic EEG FC using band-limited power (Hilbert envelope)
% using EXACT SAME FINAL fMRI bipolar seed list (after de-overlap) and SAME edge mask.

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
eegPath = fullfile(C.outDir, sprintf("EEG_segments_%s.mat", C.subjectID));
load(eegPath, "EEG");

% ---------- Load FINAL fMRI FC outputs (authoritative seedNames + pairKeepVec) ----------
fmriFCPath = fullfile(C.outDir, sprintf("FMRI_FC_%s.mat", C.subjectID));
if ~isfile(fmriFCPath)
    error("FMRI_FC file not found. Run 04_compute_fc_fmri first: %s", fmriFCPath);
end
S = load(fmriFCPath, "OUT");
OUT_FMRI = S.OUT;

% ---------- Helpers ----------
addpath(C.codeDir);
addpath(fullfile(C.codeDir, "helpers"));

% ---------- Canonical bands ----------
bands = struct();
bands.delta     = [1 3];
bands.theta     = [4 7];
bands.alpha     = [8 12];
bands.beta      = [13 30];
bands.gammalow  = [31 50];
bands.gammahigh = [51 100];

fs_ds = 250;  % envelope downsample rate

OUT = struct();

for r = 1:height(segTable)
    runID  = string(segTable.RunID(r));
    runKey = matlab.lang.makeValidName(runID);

    % ---------- Check fMRI FC exists for this run ----------
    if ~isfield(OUT_FMRI, runKey) || ~isfield(OUT_FMRI.(runKey), "seedNames")
        warning("[SKIP] %s %s: no fMRI FC seedNames found (run may have been skipped in fMRI FC).", C.subjectID, runID);
        continue;
    end

    % FINAL fMRI seed list (already de-overlapped in fMRI FC script)
    fmriSeeds_final = string(OUT_FMRI.(runKey).seedNames);
    fmriSeeds_final = upper(strrep(fmriSeeds_final, "_", "-"));
    fmriSeeds_final = canonical_bipolar_name(fmriSeeds_final);

    % exact same edge mask as fMRI FC
    if ~isfield(OUT_FMRI.(runKey), "pairKeepVec") || ~isfield(OUT_FMRI.(runKey), "nPairsKeep")
        error("[ERROR] %s %s: fMRI OUT missing pairKeepVec/nPairsKeep.", C.subjectID, runID);
    end
    pairKeepVec = OUT_FMRI.(runKey).pairKeepVec;
    nPairsKeep  = OUT_FMRI.(runKey).nPairsKeep;

    % ---------- Load EEG mono ----------
    if ~isfield(EEG, runKey)
        warning("[SKIP] %s %s: EEG struct missing this run.", C.subjectID, runID);
        continue;
    end

    eeg_mono    = EEG.(runKey).data_seg;   % [nChan x nSamp]
    fs          = EEG.(runKey).fs;
    monoLabels  = string(EEG.(runKey).labels);

    % ---------- Build bipolar EEG channels ----------
    [eeg_bip, eegBipNames] = build_bipolar_signals_from_monopolar(eeg_mono, monoLabels);

    eegBipNames = upper(strrep(eegBipNames, "_", "-"));
    eegBipNames = canonical_bipolar_name(eegBipNames);

    % ---------- Match EEG bipolars to FINAL fMRI seed list (must be exact) ----------
    tf = ismember(fmriSeeds_final, eegBipNames);

    if ~all(tf)
        missing = fmriSeeds_final(~tf);
        warning("[SKIP] %s %s: EEG missing %d/%d required fMRI seeds. (Example missing: %s)", ...
            C.subjectID, runID, nnz(~tf), numel(tf), missing(1));
        continue;
    end

    % reorder EEG to EXACT fMRI final order
    [~, idxEEG] = ismember(fmriSeeds_final, eegBipNames);
    if any(idxEEG==0)
        warning("[SKIP] %s %s: unexpected 0 indices after ismember.", C.subjectID, runID);
        continue;
    end

    eeg_bip     = eeg_bip(idxEEG, :);    % [N x nSamp]
    eegBipNames = fmriSeeds_final;       % exact labels (FMRI final order)
    nChan       = numel(eegBipNames);

    % sanity check: pairKeepVec length must match nChan choose 2
    expectedPairs = nChan*(nChan-1)/2;
    if numel(pairKeepVec) ~= expectedPairs
        error("[ERROR] %s %s: pairKeepVec length mismatch. expected=%d, got=%d", ...
            C.subjectID, runID, expectedPairs, numel(pairKeepVec));
    end

    fprintf("[OK] %s %s using FINAL fMRI seeds: N=%d | edges kept=%d\n", ...
        C.subjectID, runID, nChan, nPairsKeep);

    % ---------- Sliding windows in EEG samples ----------
    winSamp  = round(C.winSec * fs_ds);
    stepSamp = round((C.winSec - C.overlapSec) * fs_ds);

    OUT.(runKey).labels     = eegBipNames;
    OUT.(runKey).fs_raw     = fs;
    OUT.(runKey).fs_env     = fs_ds;
    OUT.(runKey).pairKeepVec = pairKeepVec;
    OUT.(runKey).nPairsKeep  = nPairsKeep;

    bandNames = fieldnames(bands);

    winStarts = [];
    winEnds   = [];
    nWin      = 0;

    for b = 1:numel(bandNames)
        bn = bandNames{b};
        fr = bands.(bn);

        % Hilbert envelope (BLP)
        env = bandpass_hilbert_envelope(eeg_bip, fs, fr(1), fr(2), fs_ds); % [nChan x nSamp_env]
        env = zscore(env, 0, 2);

        if isempty(winStarts)
            if size(env,2) < winSamp
                warning("[SKIP] %s %s: too short after ds for %ds window. band=%s", C.subjectID, runID, C.winSec, bn);
                continue;
            end
            [winStarts, winEnds] = sliding_window_indices(size(env,2), winSamp, stepSamp, 0);
            nWin = numel(winStarts);
            OUT.(runKey).winStarts = winStarts;
            OUT.(runKey).winEnds   = winEnds;
        end

        % ---------- Static EEG FC ----------
        FC_static = corr(env', "Rows","complete");   % [nChan x nChan]
        FC_static = max(min(FC_static, 0.999999), -0.999999);

        staticVec   = uppertri_vec(FC_static);
        staticVec   = staticVec(pairKeepVec);

        FC_static_z = atanh(FC_static);
        staticVec_z = uppertri_vec(FC_static_z);
        staticVec_z = staticVec_z(pairKeepVec);

        % ---------- Dynamic EEG FC ----------
        dynVec   = nan(nPairsKeep, nWin);
        dynVec_z = nan(nPairsKeep, nWin);

        for w = 1:nWin
            seg = env(:, winStarts(w):winEnds(w));
            Cw  = corr(seg', "Rows","complete");

            v = uppertri_vec(Cw);
            v = v(pairKeepVec);
            v = max(min(v, 0.999999), -0.999999);

            dynVec(:,w)   = v;
            dynVec_z(:,w) = atanh(v);
        end

        OUT.(runKey).(bn).FC_static    = FC_static;
        OUT.(runKey).(bn).FC_static_z  = FC_static_z;

        OUT.(runKey).(bn).staticVec    = staticVec;
        OUT.(runKey).(bn).staticVec_z  = staticVec_z;

        OUT.(runKey).(bn).dynVec       = dynVec;
        OUT.(runKey).(bn).dynVec_z     = dynVec_z;

        OUT.(runKey).(bn).stdDyn_z     = std(dynVec_z, 0, 2, "omitnan");
    end
end

save(fullfile(C.outDir, sprintf("EEG_FC_Hilbert_%s.mat", C.subjectID)), "OUT", "-v7.3");
fprintf("[OK] Saved %s\n", fullfile(C.outDir, sprintf("EEG_FC_Hilbert_%s.mat", C.subjectID)));
