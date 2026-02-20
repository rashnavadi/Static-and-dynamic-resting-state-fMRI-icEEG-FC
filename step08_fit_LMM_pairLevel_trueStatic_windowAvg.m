function step08_fit_LMM_pairLevel_trueStatic_windowAvg(masterDir, outDir)

% Step08_group_LMM_from_masterTables (Model A and Model B)
% Pair-level (average / trait-like) analysis
% Uses master tables built by Step07_build_master_tables to reproduce:
    % Table 2 ICC: windowAvg fMRI and windowAvg EEG per band
    % Fig3A: trueStatic EEG → trueStatic fMRI, and windowAvg EEG → windowAvg fMRI
    % Fig3B: windowAvg EEG + windowed SD → windowAvg fMRI

% INPUTS
%   Input file: master_staticFC.csv
%   outDir    : folder to save outputs (created if doesn't exist)
%
% NOTES
% - Assumes FC values in master tables are already Fisher-z.
% - ICC computed on PAIR-LEVEL dynamic values (matches manuscript "24,750 pairwise values").
% - LMM random effect: (1|subject_id)

% RUN LIKE BELOW:
% masterDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables";
% outDir    = fullfile(masterDir, "GROUP_MASTER", "group_stats_outputs_pairlevel");
% 
% step08_fit_LMM_pairLevel_trueStatic_windowAvg(masterDir, outDir)


if nargin < 2 || isempty(outDir)
    outDir = fullfile(masterDir, "group_stats_outputs");
end
if ~exist(outDir, "dir"); mkdir(outDir); end

subjDirs = dir(fullfile(masterDir, "ICE*"));
assert(~isempty(subjDirs), "No ICE* subject folders found in %s", masterDir);

T = table();
for i = 1:numel(subjDirs)
    pairFile = fullfile(subjDirs(i).folder, subjDirs(i).name, "master_staticFC.csv");
    if exist(pairFile, "file")
        Ti = readtable(pairFile);
        T  = [T; Ti]; %#ok<AGROW>
    else
        warning("Missing master_staticFC.csv for %s", subjDirs(i).name);
    end
end

assert(~isempty(T), "No pair-level data loaded.");


% ---------- Standardize band naming/order ----------
% Expected bands in master table:
bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

T.band = string(T.band);
% Just in case something sneaks in:
T.band(T.band=="gammalow")  = "low_gamma";
T.band(T.band=="gammahigh") = "high_gamma";
T.band = categorical(T.band, cellstr(bandOrder), "Ordinal", true);

% ---------- Ensure types ----------
T.subject_id = categorical(string(T.subject_id));
T.pair_id    = categorical(string(T.pair_id));
if ismember("sex", T.Properties.VariableNames)
    T.sex = categorical(string(T.sex));
else
    T.sex = categorical(repmat("NA", height(T), 1));
end

% Covariates may be missing/NaN; fitlme will drop rows with NaN in used vars.
requiredCols = ["fmri_fc_trueStatic","fmri_fc_windowAvg", ...
                "eeg_fc_trueStatic","eeg_fc_windowAvg","eeg_fc_windowed_sd", ...
                "pair_dist_mm","age"];

for c = requiredCols
    if ~ismember(c, string(T.Properties.VariableNames))
        error("Missing required column '%s' in master_staticFC.csv", c);
    end
end
T.pair_dist_mm = double(T.pair_dist_mm);
T.log_pair_dist = log(T.pair_dist_mm);
T.c_log_pair_dist = T.log_pair_dist - mean(T.log_pair_dist, 'omitnan');


% ---------- 1) TABLE 2 ICC ----------
Ticc = compute_table2_icc(T, bandOrder);
writetable(Ticc, fullfile(outDir, "Table2_ICC_windowAvg.csv"));
save(fullfile(outDir, "Table2_ICC_windowAvg.mat"), "Ticc", "-v7.3");

% ---------- 2) FIG 3A LMMs (static, per band) ----------
Tout_trueStatic  = run_lmm_per_band(T, bandOrder, ...
    "fmri_fc_trueStatic", "eeg_fc_trueStatic", ...
    "fmri_fc_trueStatic ~ eeg_fc_trueStatic + eeg_fc_windowed_sd + c_log_pair_dist + age + sex + (1 + eeg_fc_trueStatic | subject_id) + (1 | pair_id)");

Tout_windowAvg = run_lmm_per_band(T, bandOrder, ...
    "fmri_fc_windowAvg", "eeg_fc_windowAvg", ...
    "fmri_fc_windowAvg ~ eeg_fc_windowAvg + eeg_fc_windowed_sd + c_log_pair_dist + age + sex + (1 + eeg_fc_windowAvg | subject_id) + (1 | pair_id)");

% FDR across 6 bands within each family
Tout_trueStatic.p_fdr  = fdr_bh(Tout_trueStatic.p_unc);
Tout_windowAvg.p_fdr   = fdr_bh(Tout_windowAvg.p_unc);

writetable(Tout_trueStatic, fullfile(outDir, "Fig3A_LMM_trueStatic.csv"));
writetable(Tout_windowAvg,  fullfile(outDir, "Fig3A_LMM_windowAvg.csv"));

save(fullfile(outDir, "Fig3A_LMM_trueStatic.mat"),  "Tout_trueStatic", "-v7.3");
save(fullfile(outDir, "Fig3A_LMM_windowAvg.mat"),   "Tout_windowAvg",  "-v7.3");

% ---------- 3) FIG 3B LMMs (dynamic + SD, per band) ----------
[Tout_dynTerm, Tout_sdTerm] = run_lmm_dynamic_with_sd(T, bandOrder);

Tout_dynTerm.p_fdr = fdr_bh(Tout_dynTerm.p_unc);
Tout_sdTerm.p_fdr  = fdr_bh(Tout_sdTerm.p_unc);

writetable(Tout_dynTerm, fullfile(outDir, "Fig3B_LMM_windowAvgTerm.csv"));
writetable(Tout_sdTerm,  fullfile(outDir, "Fig3B_LMM_windowedSDterm.csv"));
save(fullfile(outDir, "Fig3B_LMM_windowAvgTerm.mat"),   "Tout_dynTerm", "-v7.3");
save(fullfile(outDir, "Fig3B_LMM_windowedSDterm.mat"),  "Tout_sdTerm",  "-v7.3");

fprintf("\n[OK] Group stats complete.\nSaved outputs to: %s\n", outDir);

end

% ===================== Helpers =====================

function Ticc = compute_table2_icc(T, bandOrder)
% ICC computed on PAIR-LEVEL window-averaged values (mean of windowed FC), grouped by subject_id.

n = 1 + numel(bandOrder);  % fMRI + EEG bands
Ticc = table('Size',[n 3], ...
    'VariableTypes', {'string','double','double'}, ...
    'VariableNames', {'Modality','ICC','SD'});

k = 1;

% fMRI ICC on dynamic (use all bands, not only delta)
[iccVal, sdVal] = icc_oneway(T.fmri_fc_windowAvg, T.subject_id);

Ticc.Modality(k) = "fMRI";
Ticc.ICC(k) = iccVal;
Ticc.SD(k)  = sdVal;
k = k + 1;

% EEG ICC per band on dynamic
for b = bandOrder
    Tb = T(string(T.band)==b, :);
    [iccVal, sdVal] = icc_oneway(Tb.eeg_fc_windowAvg, Tb.subject_id);
    Ticc.Modality(k) = string(b);
    Ticc.ICC(k) = iccVal;
    Ticc.SD(k)  = sdVal;
    k = k + 1;
end
end


function [iccVal, sdPlaceholder] = icc_oneway(y, subj)
% ICC(1) via one-way random effects ANOVA components.
% sdPlaceholder is not true ICC SD (manuscript SD may have been bootstrap).
% If you have old bootstrap code, replace this placeholder.
y = double(y);
subj = categorical(subj);

ok = ~isnan(y) & ~ismissing(subj);
y = y(ok); subj = subj(ok);

if numel(unique(subj)) < 2
    iccVal = NaN;
    sdPlaceholder = NaN;
    return;
end

[~, tbl] = anova1(y, subj, "off");
MSb = tbl{2,4};
MSw = tbl{3,4};

counts = countcats(subj);
kbar = mean(counts);

iccVal = (MSb - MSw) / (MSb + (kbar - 1) * MSw);

% Placeholder (set to NaN to avoid implying correctness)
sdPlaceholder = NaN;
end

function Tout = run_lmm_per_band(T, bandOrder, yCol, xCol, formula)
% Fit per band and extract fixed effect for xCol

rows = {};
for b = bandOrder
    Tb = T(string(T.band)==b, :);

    % Remove rows with NaNs in required columns
    useCols = string([yCol, xCol, "age", "sex", "subject_id"]);
    Tb = rmmissing(Tb, "DataVariables", useCols);

    if height(Tb) < 10 || numel(categories(Tb.subject_id)) < 2
        rows(end+1,:) = {b, NaN, NaN, NaN, NaN, height(Tb)};
        continue;
    end

    lme = fitlme(Tb, formula, "FitMethod","REML");
    r = lmm_extract(lme, xCol);

    rows(end+1,:) = {b, r.beta, r.p, r.ci_low, r.ci_high, height(Tb)};
end

Tout = cell2table(rows, "VariableNames", ...
    ["band","beta","p_unc","ci_low","ci_high","nObs"]);
end

function [Tout_dynTerm, Tout_sdTerm] = run_lmm_dynamic_with_sd(T, bandOrder)
% Per band: fmri_fc_windowAvg ~ eeg_fc_windowAvg + eeg_fc_windowed_sd + age + sex + (1|subject_id)
% Extract both terms separately

rowsDyn = {};
rowsSD  = {};
formula = "fmri_fc_windowAvg ~ eeg_fc_windowAvg + eeg_fc_windowed_sd + c_log_pair_dist + age + sex + (1|subject_id) + (1|pair_id)";

for b = bandOrder
    Tb = T(string(T.band)==b, :);

    useCols = ["fmri_fc_windowAvg","eeg_fc_windowAvg","eeg_fc_windowed_sd","age","sex","subject_id"];

    Tb = rmmissing(Tb, "DataVariables", useCols);

    if height(Tb) < 10 || numel(categories(Tb.subject_id)) < 2
        rowsDyn(end+1,:) = {b, NaN, NaN, NaN, NaN, height(Tb)};
        rowsSD(end+1,:)  = {b, NaN, NaN, NaN, NaN, height(Tb)};
        continue;
    end

    lme = fitlme(Tb, formula, "FitMethod","REML");

    rd = lmm_extract(lme, "eeg_fc_windowAvg");
    rs = lmm_extract(lme, "eeg_fc_windowed_sd");

    rowsDyn(end+1,:) = {b, rd.beta, rd.p, rd.ci_low, rd.ci_high, height(Tb)};
    rowsSD(end+1,:)  = {b, rs.beta, rs.p, rs.ci_low, rs.ci_high, height(Tb)};
end

Tout_dynTerm = cell2table(rowsDyn, "VariableNames", ...
    ["band","beta","p_unc","ci_low","ci_high","nObs"]);

Tout_sdTerm  = cell2table(rowsSD, "VariableNames", ...
    ["band","beta","p_unc","ci_low","ci_high","nObs"]);
end

function r = lmm_extract(lme, termName)
coef = lme.Coefficients;
idx = strcmp(coef.Name, termName);

if ~any(idx)
    error("Term '%s' not found in model coefficients.", termName);
end

beta = coef.Estimate(idx);
p    = coef.pValue(idx);
ciLo = coef.Lower(idx);
ciHi = coef.Upper(idx);

r.beta = beta;
r.p = p;
r.ci_low = ciLo;
r.ci_high = ciHi;
end

function p_adj = fdr_bh(p)
% Benjamini–Hochberg adjusted p-values
p = double(p(:));
p(isnan(p)) = 1; % treat NaNs as non-significant

[ps, idx] = sort(p);
m = numel(p);
q = ps .* m ./ (1:m)';

% enforce monotonicity
q = flipud(cummin(flipud(q)));
q(q > 1) = 1;

p_adj = nan(size(p));
p_adj(idx) = q;
end
