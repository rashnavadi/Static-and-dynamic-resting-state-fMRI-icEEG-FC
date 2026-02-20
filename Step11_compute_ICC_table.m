% Step09_compute_ICC_table
% Computes ICC (mean ± SD via subject bootstrap) for:
%   1) dynamic fMRI FC (window-level)
%   2) dynamic EEG FC (window-level) per frequency band
%
% Inputs
%   masterOutDir : folder that contains subject subfolders written by Step07
%                  e.g., masterOutDir/ICE013/master_dynamicFC.csv
%   outDir       : where to save ICC table CSV
%   nBoot        : number of bootstrap resamples (default: 500)
%
% Output
%   ICC_table_window_level.csv
% output saved in this dir: GROUP_MASTER/group_stats_outputs/ICC_table_window_level.csv

% RUN example:
% masterOutDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables";
% outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs";
% 
% Step11_compute_ICC_table(masterOutDir, outDir, 500);

function Step11_compute_ICC_table(masterOutDir, outDir, nBoot)
% Computes ICC (mean ± SD via subject bootstrap) for:
%   1) dynamic fMRI FC (window-level): fmri_fc_windowed
%   2) dynamic EEG FC (window-level) per band: eeg_fc_windowed
%
% Reads: masterOutDir/ICE###/master_dynamicFC.csv
%
% Output: ICC_table_window_level.csv

if nargin < 3 || isempty(nBoot), nBoot = 500; end
if ~exist(outDir, "dir"), mkdir(outDir); end

bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

% -------------------------
% 1) Load & concatenate all subjects
% -------------------------
subDirs = dir(fullfile(masterOutDir, "ICE*"));
subDirs = subDirs([subDirs.isdir]);

Tall = table();
nLoaded = 0;

for i = 1:numel(subDirs)
    csvPath = fullfile(subDirs(i).folder, subDirs(i).name, "master_dynamicFC.csv");
    if ~exist(csvPath, "file"), continue; end

    T = readtable(csvPath);
    v = string(T.Properties.VariableNames);

    req = ["subject_id","band","fmri_fc_windowed","eeg_fc_windowed","run_id","pair_id","win_idx"];
    if any(~ismember(req, v))
        warning("Missing required columns in %s (skipping). Found: %s", csvPath, strjoin(v,", "));
        continue;
    end

    % Normalize types
    T.subject_id = string(T.subject_id);
    T.run_id     = string(T.run_id);
    T.pair_id    = string(T.pair_id);
    T.band       = string(T.band);

    Tall = [Tall; T]; %#ok<AGROW>
    nLoaded = nLoaded + 1;
end

if nLoaded == 0 || isempty(Tall)
    error("No usable master_dynamicFC.csv files found in %s", masterOutDir);
end

subjects = unique(Tall.subject_id);
nSubj = numel(subjects);
fprintf("[OK] Loaded %d subjects (%d rows).\n", nLoaded, height(Tall));

% -------------------------
% 2) fMRI ICC (DEDUPLICATE across bands)
% -------------------------
% fMRI does not depend on EEG band. In this table, fMRI is repeated for each band.
% Deduplicate by unique (subject_id, run_id, pair_id, win_idx).
TfMRI = Tall(:, ["subject_id","run_id","pair_id","win_idx","fmri_fc_windowed"]);
TfMRI = unique(TfMRI, "rows");

icc_fmri = compute_icc_from_lme(TfMRI, "fmri_fc_windowed");

icc_fmri_boot = nan(nBoot,1);
for bb = 1:nBoot
    bootSubj = subjects(randi(nSubj, nSubj, 1));
    Tboot = TfMRI(ismember(string(TfMRI.subject_id), string(bootSubj)), :);
    icc_fmri_boot(bb) = compute_icc_from_lme(Tboot, "fmri_fc_windowed");
end
icc_fmri_sd = std(icc_fmri_boot, "omitnan");

% -------------------------
% 3) EEG ICC per band
% -------------------------
icc_eeg    = nan(numel(bandOrder),1);
icc_eeg_sd = nan(numel(bandOrder),1);

for k = 1:numel(bandOrder)
    bn = bandOrder(k);

    Teeg = Tall(Tall.band == bn, ["subject_id","eeg_fc_windowed"]);
    icc_eeg(k) = compute_icc_from_lme(Teeg, "eeg_fc_windowed");

    icc_boot = nan(nBoot,1);
    for bb = 1:nBoot
        bootSubj = subjects(randi(nSubj, nSubj, 1));
        Tboot = Teeg(ismember(string(Teeg.subject_id), string(bootSubj)), :);
        icc_boot(bb) = compute_icc_from_lme(Tboot, "eeg_fc_windowed");
    end
    icc_eeg_sd(k) = std(icc_boot, "omitnan");
end

% -------------------------
% 4) Save output
% -------------------------
RowNames = ["fMRI"; bandOrder(:)];
ICC = [icc_fmri; icc_eeg];
SD  = [icc_fmri_sd; icc_eeg_sd];

Tout = table(RowNames, ICC, SD, 'VariableNames', {'Modality_or_Band','ICC','SD'});

csvOut = fullfile(outDir, "ICC_table_window_level.csv");
writetable(Tout, csvOut);

fprintf("[SAVED] %s\n", csvOut);

end

% ===== helper: compute ICC from random-intercept LME =====
function icc = compute_icc_from_lme(Tin, yVar)

Tin = Tin(~isnan(Tin.(yVar)), :);
if height(Tin) < 10
    icc = NaN; return;
end

Tin.subject_id = categorical(string(Tin.subject_id));

lme = fitlme(Tin, sprintf("%s ~ 1 + (1|subject_id)", yVar), "FitMethod","REML");

% subject random-intercept variance
varSubj = NaN;
try
    varSubj = lme.CovarianceParameters.Estimate(1);
catch
    cp = covarianceParameters(lme);
    if istable(cp) && ismember("Estimate", string(cp.Properties.VariableNames))
        varSubj = cp.Estimate(1);
    elseif isnumeric(cp)
        varSubj = cp(1);
    else
        varSubj = cp{1};
        if isstruct(varSubj) && isfield(varSubj,"Estimate")
            varSubj = varSubj.Estimate;
        end
    end
end

varRes = lme.MSE;
icc = varSubj / (varSubj + varRes);

end

