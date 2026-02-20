function SupTable_compare_beta2_Model1_vs_Model2(outDir)
% SupTable_compare_beta2_Model1_vs_Model2
% Creates a supplementary CSV comparing beta2 (DEV term) between:
%   Model 1: random intercepts
%   Model 2: random slope for eeg_dyn_dev (sensitivity)
%
% OUTPUTS:
%   SupTable_beta2_Model1_vs_Model2.csv
%
% RUN:
% outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window";
% SupTable_compare_beta2_Model1_vs_Model2(outDir);

f1 = fullfile(outDir, "LMM_window_DEV_randomIntercepts_bonfAcrossBands.csv");
f2 = fullfile(outDir, "LMM_window_DEV_randomSlopeTry_bonfAcrossBands.csv");
assert(exist(f1,"file")==2, "Missing file: %s", f1);
assert(exist(f2,"file")==2, "Missing file: %s", f2);

T1 = readtable(f1);
T2 = readtable(f2);

bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

% enforce band order + sort
T1.band = categorical(string(T1.band), cellstr(bandOrder), "Ordinal", true);
T2.band = categorical(string(T2.band), cellstr(bandOrder), "Ordinal", true);
T1 = sortrows(T1,"band");
T2 = sortrows(T2,"band");

% align by band safely
[commonBands, ia, ib] = intersect(T1.band, T2.band, "stable");
T1 = T1(ia,:);
T2 = T2(ib,:);

% extract beta2 info
beta2_m1 = T1.beta_dev;
ciL_m1   = T1.ci_low_dev;
ciH_m1   = T1.ci_high_dev;
p_m1     = T1.p_dev;
pbonf_m1 = pickcol(T1, "p_dev_bonf", NaN(size(p_m1)));

beta2_m2 = T2.beta_dev;
ciL_m2   = T2.ci_low_dev;
ciH_m2   = T2.ci_high_dev;
p_m2     = T2.p_dev;
pbonf_m2 = pickcol(T2, "p_dev_bonf", NaN(size(p_m2)));

% deltas (Model2 - Model1)
d_beta = beta2_m2 - beta2_m1;

% convergence indicator (Model2 is NaN when fit failed in your script)
model2_ok = ~isnan(beta2_m2) & ~isnan(p_m2);

Tout = table();
Tout.band = commonBands;

Tout.beta2_M1      = beta2_m1;
Tout.ciLow_M1      = ciL_m1;
Tout.ciHigh_M1     = ciH_m1;
Tout.p_M1          = p_m1;
Tout.pBonf_M1      = pbonf_m1;

Tout.beta2_M2      = beta2_m2;
Tout.ciLow_M2      = ciL_m2;
Tout.ciHigh_M2     = ciH_m2;
Tout.p_M2          = p_m2;
Tout.pBonf_M2      = pbonf_m2;

Tout.delta_beta2   = d_beta;
Tout.model2_ok     = model2_ok;

% also keep sample size if present
Tout.nObs  = pickcol(T1,"nObs",  NaN(height(T1),1));
Tout.nSubj = pickcol(T1,"nSubj", NaN(height(T1),1));
Tout.nPairs= pickcol(T1,"nPairs",NaN(height(T1),1));

outCSV = fullfile(outDir, "SupTable_beta2_Model1_vs_Model2.csv");
writetable(Tout, outCSV);

fprintf("[OK] Wrote: %s\n", outCSV);

end

% ---------- helper: safe column grab ----------
function v = pickcol(T, name, defaultVal)
if ismember(name, T.Properties.VariableNames)
    v = T.(name);
else
    v = defaultVal;
end
end
