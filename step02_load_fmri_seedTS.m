%% written by Tahereh Rashnavadi
% Jan 2026



%% 02_load_fmri_seedTS_ICE033.m
% Loads the bipolar fMRI seed time series produced by:
%   extract_bipolar_fMRI_seedTS_excludeIED(...)
% Cleans invalid seeds, canonicalizes names, z-scores across time.

% clear; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

% ---------- Load config + segTable ----------
if ~exist("C","var") || ~isfield(C,"subjectID") || strlength(C.subjectID)==0
    error("C.subjectID must be set by the calling script");
end

cfgPath = fullfile(C.projectRoot, "FC_rebuild", C.subjectID, sprintf("config_%s.mat", C.subjectID));
load(cfgPath, "C");

% ---------- Helpers ----------
addpath(C.codeDir);
addpath(fullfile(C.codeDir, "helpers"));

FMRI = struct();

for r = 1:height(segTable)
    runID = string(segTable.RunID(r));
    runKey = matlab.lang.makeValidName(runID);

    % File created by extraction stage
    bipFile = fullfile(C.fmriTSDir, runID, ...
        sprintf(C.bipTS_pattern, C.subjectID, runID));

    if ~isfile(bipFile)
        error("Missing fMRI bipolar TS file: %s", bipFile);
    end

    S = load(bipFile);
    seedTS_bipolar   = S.seedTS_bipolar;     % [T x Nseed]
    seedNames_bipolar= string(S.seedNames_bipolar);

    % ---------- Normalize names ----------
    seedNames_bipolar = upper(strrep(seedNames_bipolar, "_", "-"));
    seedNames_bipolar = canonical_bipolar_name(seedNames_bipolar);

    % ---------- Remove invalid seeds ----------
    bad = all(seedTS_bipolar==0,1) | all(~isfinite(seedTS_bipolar),1) | (std(seedTS_bipolar,0,1,'omitnan')==0);
    if any(bad)
        fprintf("[INFO] %s %s removing %d invalid fMRI seeds (zero/NaN/flat)\n", C.subjectID, runID, nnz(bad));
        seedTS_bipolar(:,bad) = [];
        seedNames_bipolar(bad)= [];
        if isfield(S,"D") && ~isempty(S.D)
            D = S.D;
            D(bad,:) = []; D(:,bad) = [];
        else
            D = [];
        end
    else
        if isfield(S,"D"), D = S.D; else, D = []; end
    end

    % ---------- Z-score across time ----------
    seedTS_bipolar = zscore(seedTS_bipolar, 0, 1);

    % ---------- Store ----------
    FMRI.(runKey).seedTS    = seedTS_bipolar;
    FMRI.(runKey).seedNames = seedNames_bipolar;
    FMRI.(runKey).runID     = runID;   % keep original run label (optional but helpful)

    if isfield(S,"voxelSize"), FMRI.(runKey).voxelSize = S.voxelSize; else, FMRI.(runKey).voxelSize = []; end
    FMRI.(runKey).D = D;

    if isfield(S,"excludedSeeds"),  FMRI.(runKey).excludedSeeds  = string(S.excludedSeeds);  end
    if isfield(S,"excludedReason"), FMRI.(runKey).excludedReason = string(S.excludedReason); end

    fprintf("[OK] Loaded fMRI seeds: %s %s (T=%d, N=%d)\n", C.subjectID, runID, size(seedTS_bipolar,1), size(seedTS_bipolar,2));
end

save(fullfile(C.outDir, sprintf("FMRI_seedTS_%s.mat", C.subjectID)), "FMRI", "-v7.3");
fprintf("[OK] Saved %s\n", fullfile(C.outDir, sprintf("FMRI_seedTS_%s.mat", C.subjectID)));
