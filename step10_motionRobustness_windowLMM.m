%% STEP10 Motion robustness: low-motion vs moderate-motion runs (FEAT rel.rms)
% Purpose:
%   Show EEG–fMRI coupling is preserved even in low-motion runs.
% Strategy:
%   1) Compute run-level motion metric from FEAT: mc/prefiltered_func_data_mcf_rel.rms
%   2) Join run-level motion into window-level LMM table (Tall)
%   3) Stratify runs into low vs moderate motion
%   4) Fit the same LMM separately in each subgroup (band-resolved)
%   5) Try lowThresh=0.2 mm first; if low-motion group too small => fallback to 0.5 mm
%
% How to run:
%   step10_motionRobustness_windowLMM;

function step10_motionRobustness_windowLMM()

%% -------------------- USER PATHS --------------------
partsDir   = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/parts_window";
preprocRoot= "/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis"; % contains ICE*/.../*.feat
outDir     = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_motion";
if ~exist(outDir,"dir"); mkdir(outDir); end

% Optional caching (recommended)
cacheMotionMat = fullfile(outDir, "TrunMotion_feat_relRMS.mat");
rebuildMotionTable = true; % set true if FEAT outputs changed

%% -------------------- THRESHOLDS --------------------
% - exclude anything above 1.5 mm (should already be handled at window level)
% - split remaining runs into low vs moderate based on run-level max motion
tryLowThresh = 0.2;  % try this first
fallbackThresh = 0.5; % if low group too small
excludeCeil = 1.5;   % ceiling for "moderate" group

% Minimum data requirements (same spirit as your main/spatial)
minObs   = 2000;
minSubj  = 10;   % Pierre’s suggestion for "enough subjects"
minRuns  = 12;   % practical safeguard
minPairs = 50;

%% -------------------- MODEL --------------------
bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

formula = ...
 "fmri_fc_windowed ~ eeg_fc_trueStatic + delta_eeg_fc + eeg_fc_windowed_sd + c_log_pair_dist + age + sex" + ...
 " + (1|subject_id) + (1|subject_id:run_id) + (1|pair_id)";

%% ====================================================
%  A) Build / Load run-level motion table (TrunMotion)
% =====================================================
if ~rebuildMotionTable && exist(cacheMotionMat,"file")
    load(cacheMotionMat, "TrunMotion");
    fprintf("[OK] Loaded cached motion table: %s\n", cacheMotionMat);
else
    TrunMotion = build_motion_table_feat_relRMS(preprocRoot);
    save(cacheMotionMat, "TrunMotion", "-v7.3");
    fprintf("[OK] Built + cached motion table: %s\n", cacheMotionMat);
end

assert(all(ismember(["subject_id","run_id","maxDisp_run"], TrunMotion.Properties.VariableNames)), ...
    "TrunMotion must contain subject_id, run_id, maxDisp_run.");

%% ====================================================
%  B) Load + concat window-level tables (Tall)
% =====================================================
files = dir(fullfile(partsDir, "ICE*_forLMM_window.mat"));
assert(~isempty(files), "No *_forLMM_window.mat found in %s", partsDir);

Tall = table();
for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name), "Tsub");
    Tall = [Tall; S.Tsub]; %#ok<AGROW>
end

% Basic type fixes (match your style)
Tall.subject_id = categorical(string(Tall.subject_id));
Tall.run_id     = categorical(string(Tall.run_id));
Tall.pair_id    = string(Tall.pair_id);
Tall.sex        = categorical(string(Tall.sex));
Tall.band       = categorical(string(Tall.band), cellstr(bandOrder), "Ordinal", true);

disp("Tall subject examples:");
disp(unique(string(Tall.subject_id(1:20)))');

disp("TrunMotion subject examples:");
disp(unique(string(TrunMotion.subject_id(1:20)))');
%% ====================================================
%  C) Join motion into Tall
% =====================================================
% --- CREATE JOIN KEY IN BOTH TABLES (needed even when loading cached TrunMotion) ---
Tall.run_join = normalize_run_id(Tall.run_id);

if ~ismember("run_join", TrunMotion.Properties.VariableNames)
    TrunMotion.run_join = normalize_run_id(TrunMotion.run_id);
end

fprintf("[QC] Example Tall run_id -> run_join:\n");
disp(Tall(1:10, {'run_id','run_join'}));

fprintf("[QC] Example TrunMotion run_id -> run_join:\n");
disp(TrunMotion(1:10, {'run_id','run_join'}));

% Avoid variable name collision in outerjoin (Tall already has run_id)
if ismember("run_id", TrunMotion.Properties.VariableNames)
    TrunMotion = renamevars(TrunMotion, "run_id", "run_id_feat");
end

Tall = outerjoin(Tall, TrunMotion, ...
    "Keys", ["subject_id","run_join"], ...
    "MergeKeys", true, ...
    "Type","left");


if ~ismember("maxDisp_run", Tall.Properties.VariableNames)
    error("Join failed: maxDisp_run not added. Check subject_id/run_id consistency.");
end

% ---------------------------
% Drop rows with missing motion (runs not found)
% ---------------------------
missMotion = isnan(Tall.maxDisp_run);

nBefore = height(Tall);
nMiss   = sum(missMotion);

fprintf("[QC] Missing maxDisp_run after join: %d of %d rows (%.2f%%)\n", ...
    nMiss, nBefore, 100*nMiss/nBefore);

if nMiss > 0
    fprintf("[WARN] Dropping %d rows missing maxDisp_run (unmatched subject/run)\n", nMiss);
    Tall(missMotion,:) = [];
end

fprintf("[QC] Rows before drop: %d\n", nBefore);
fprintf("[QC] Rows remaining: %d\n", height(Tall));


%% ====================================================
%  D) Define motion group with 0.2 threshold first
% =====================================================
Tall = label_motion_groups(Tall, tryLowThresh, excludeCeil);

% QC summary (overall)
fprintf("\n========== QC: Motion group counts (initial threshold %.1f mm) ==========\n", tryLowThresh);
print_motion_qc(Tall, bandOrder);

% Decide if threshold is acceptable (low group big enough)
ok = motion_threshold_ok(Tall, bandOrder, minSubj, minRuns);

if ~ok
    fprintf("\n[INFO] Low-motion group too small at %.1f mm. Falling back to %.1f mm.\n", tryLowThresh, fallbackThresh);
    Tall = label_motion_groups(Tall, fallbackThresh, excludeCeil);

    fprintf("\n========== QC: Motion group counts (fallback threshold %.1f mm) ==========\n", fallbackThresh);
    print_motion_qc(Tall, bandOrder);

    usedThresh = fallbackThresh;
else
    usedThresh = tryLowThresh;
end

% Keep only low/moderate
Tall = Tall(Tall.motion_group=="low" | Tall.motion_group=="moderate", :);
Tall.motion_group = removecats(Tall.motion_group);

%% ====================================================
%  E) Ensure distance covariate exists
% =====================================================
if ~ismember("c_log_pair_dist", Tall.Properties.VariableNames)
    assert(ismember("pair_dist_mm", Tall.Properties.VariableNames), "No distance column found.");
    Tall.pair_dist_mm = double(Tall.pair_dist_mm);
    Tall.log_pair_dist = log(Tall.pair_dist_mm);
    Tall.c_log_pair_dist = Tall.log_pair_dist - mean(Tall.log_pair_dist, 'omitnan');
end

%% ====================================================
%  F) Run LMM per motion group × band
% =====================================================
Results = table();
groups = categories(Tall.motion_group);

fprintf("\n========== Fitting LMMs (threshold %.1f mm) ==========\n", usedThresh);

for gi = 1:numel(groups)
    grp = groups{gi};
    Tg0 = Tall(Tall.motion_group==grp, :);

    for b = bandOrder
        Tb = Tg0(string(Tg0.band)==string(b), :);

        useCols = ["fmri_fc_windowed","eeg_fc_trueStatic","delta_eeg_fc","eeg_fc_windowed_sd", ...
                   "c_log_pair_dist","age","sex","subject_id","run_id","pair_id"];
        Tb = rmmissing(Tb, "DataVariables", useCols);

        nObs  = height(Tb);
        nSubj = numel(unique(string(Tb.subject_id)));
        nRuns = numel(unique(string(Tb.subject_id) + "_" + string(Tb.run_id)));
        nPairs= numel(unique(string(Tb.pair_id)));

        fprintf("[QC] %s | %s: nObs=%d nSubj=%d nRuns=%d nPairs=%d\n", grp, string(b), nObs, nSubj, nRuns, nPairs);

        if nObs < minObs || nSubj < minSubj || nRuns < minRuns || nPairs < minPairs
            fprintf("     -> [SKIP] below thresholds\n");
            continue;
        end

        try
            lme = fitlme(Tb, formula, "FitMethod","REML");

            rDev = lmm_extract(lme, "delta_eeg_fc");
            rSta = lmm_extract(lme, "eeg_fc_trueStatic");

            Results = [Results; table( ...
                string(grp), string(b), ...
                rDev.beta, rDev.p, rDev.ci_low, rDev.ci_high, ...
                rSta.beta, rSta.p, rSta.ci_low, rSta.ci_high, ...
                nObs, nSubj, nRuns, nPairs, usedThresh, ...
                'VariableNames', ["motion_group","band", ...
                                  "beta_dev","p_dev","ci_low_dev","ci_high_dev", ...
                                  "beta_static","p_static","ci_low_static","ci_high_static", ...
                                  "nObs","nSubj","nRuns","nPairs","lowThresh_mm"])]; %#ok<AGROW>
        catch ME
            fprintf("[WARN] LMM failed for %s | %s: %s\n", grp, string(b), ME.message);
        end
    end
end

% Save results
csvOut = fullfile(outDir, "MotionRobustness_LMM_byGroup_thr" + strrep(num2str(usedThresh),".","p") + ".csv");
matOut = fullfile(outDir, "MotionRobustness_LMM_byGroup_thr" + strrep(num2str(usedThresh),".","p") + ".mat");

writetable(Results, csvOut);
save(matOut, "Results", "-v7.3");

fprintf("\n[OK] Saved motion robustness results:\n%s\n", csvOut);

% (Optional) quick plot helper: effect sizes by band & group
% make_motion_robustness_plot(Results, outDir);

end % main function


%% ===================== HELPERS =====================

function TrunMotion = build_motion_table_feat_relRMS(preprocRoot)
% Build run-level motion metrics from FEAT:
%   *.feat/mc/prefiltered_func_data_mcf_rel.rms

featDirs = dir(fullfile(preprocRoot, "ICE*", "4_MRI", "2_Functionals", "2_PreProcessing", "*.feat"));
assert(~isempty(featDirs), "No *.feat directories found under %s", preprocRoot);

S = struct('subject_id',{},'run_id',{},'maxDisp_run',{},'meanDisp_run',{},'p95Disp_run',{},'nFrames',{});

for i = 1:numel(featDirs)
    featPath = fullfile(featDirs(i).folder, featDirs(i).name);

    % subject_id from path
    parts = split(string(featPath), filesep);

    % Find the first token that matches ICE + digits (e.g., ICE013)
    idx = find(~cellfun('isempty', regexp(cellstr(parts), '^ICE\d+$', 'once')), 1, 'first');

    if isempty(idx)
        fprintf("[SKIP] Could not parse subject_id from path: %s\n", featPath);
        continue;
    end

    subj = parts(idx);


    % run_id from folder name like "Run1a.feat"
    run_id = erase(string(featDirs(i).name), ".feat");
    % Skip excluded FEAT folders (they were not used downstream)
    if contains(run_id, "Exclude", "IgnoreCase", true)
        continue;
    end

    rmsPath = fullfile(featPath, "mc", "prefiltered_func_data_mcf_rel.rms");
    if ~exist(rmsPath,"file")
        cand = dir(fullfile(featPath, "mc", "*rel*.rms"));
        if isempty(cand)
            fprintf("[SKIP] Missing rel.rms: %s\n", rmsPath);
            continue;
        end
        rmsPath = fullfile(cand(1).folder, cand(1).name);
    end

    disp_mm = readmatrix(rmsPath, "FileType", "text");
    disp_mm = disp_mm(:);
    disp_mm = disp_mm(~isnan(disp_mm));

    if isempty(disp_mm)
        fprintf("[SKIP] Empty rel.rms: %s\n", rmsPath);
        continue;
    end

    S(end+1).subject_id  = string(subj);
    S(end).run_id        = string(run_id);
    S(end).maxDisp_run   = max(disp_mm);
    S(end).meanDisp_run  = mean(disp_mm);
    S(end).p95Disp_run   = prctile(disp_mm,95);
    S(end).nFrames       = numel(disp_mm);
end

TrunMotion = struct2table(S);
TrunMotion.subject_id = categorical(TrunMotion.subject_id);
TrunMotion.run_id     = categorical(TrunMotion.run_id);

end


function Tall = label_motion_groups(Tall, lowThresh, excludeCeil)
% Label runs into low/moderate/excluded using run-level maxDisp_run (mm)
Tall.motion_group = strings(height(Tall),1);

Tall.motion_group(Tall.maxDisp_run < lowThresh) = "low";
Tall.motion_group(Tall.maxDisp_run >= lowThresh & Tall.maxDisp_run <= excludeCeil) = "moderate";
Tall.motion_group(Tall.maxDisp_run > excludeCeil) = "excluded";

Tall.motion_group = categorical(Tall.motion_group);
Tall.motion_group = removecats(Tall.motion_group);

end


function print_motion_qc(Tall, bandOrder)
% Quick QC counts per group and per band
groups = categories(Tall.motion_group);

for gi = 1:numel(groups)
    grp = groups{gi};
    Tg = Tall(Tall.motion_group==grp, :);

    nObs_all  = height(Tg);
    nSubj_all = numel(unique(string(Tg.subject_id)));
    nRuns_all = numel(unique(string(Tg.subject_id) + "_" + string(Tg.run_id)));

    fprintf("Group=%s: nObs=%d | nSubj=%d | nRuns=%d\n", grp, nObs_all, nSubj_all, nRuns_all);

    for b = bandOrder
        Tb = Tg(string(Tg.band)==string(b), :);
        fprintf("   %s: nObs=%d nSubj=%d nRuns=%d\n", string(b), height(Tb), ...
            numel(unique(string(Tb.subject_id))), ...
            numel(unique(string(Tb.subject_id) + "_" + string(Tb.run_id))));
    end
end
end


function ok = motion_threshold_ok(Tall, bandOrder, minSubj, minRuns)
% Decide if low-motion group is big enough overall and across bands
ok = true;

Tlow = Tall(Tall.motion_group=="low", :);
if isempty(Tlow)
    ok = false; return;
end

nSubj = numel(unique(string(Tlow.subject_id)));
nRuns = numel(unique(string(Tlow.subject_id) + "_" + string(Tlow.run_id)));

if nSubj < minSubj || nRuns < minRuns
    ok = false; return;
end

% Also require that low group has at least minSubj in a majority of bands
goodBands = 0;
for b = bandOrder
    Tb = Tlow(Tlow.band==categorical(string(b)), :);
    nb = numel(unique(string(Tb.subject_id)));
    if nb >= minSubj
        goodBands = goodBands + 1;
    end
end

if goodBands < ceil(numel(bandOrder)/2)
    ok = false;
end

end


function r = lmm_extract(lme, termName)
coef = lme.Coefficients;
idx = strcmp(coef.Name, termName);
if ~any(idx)
    error("Term '%s' not found in model coefficients.", termName);
end
r.beta    = coef.Estimate(idx);
r.p       = coef.pValue(idx);
r.ci_low  = coef.Lower(idx);
r.ci_high = coef.Upper(idx);
end


function run_join = normalize_run_id(run_id)
% Convert run IDs to a canonical "Run<number>" label
% Examples:
%   "Run1a" -> "Run1"
%   "Run1b" -> "Run1"
%   "Run01" -> "Run1"
%   "run2"  -> "Run2"

r = string(run_id);
r = upper(r);

% Extract "RUN" + digits anywhere in the string
tok = regexp(r, "RUN\s*0*(\d+)", "tokens", "once");

run_join = strings(size(r));
for i = 1:numel(r)
    if isempty(tok{i})
        run_join(i) = "";   % will become missing
    else
        run_join(i) = "Run" + tok{i}{1};
    end
end

run_join = categorical(run_join);
end

