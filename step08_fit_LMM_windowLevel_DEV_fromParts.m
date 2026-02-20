function step08_fit_LMM_windowLevel_DEV_fromParts(partsDir, outDir)
% Fits WINDOW-level LMM per band with BH-FDR correction across 6 bands.
%
% DEFAULT (single) model per band (Pierre basis):
%   fmri_fc_windowed ~ eeg_fc_trueStatic + delta_eeg_fc + eeg_fc_windowed_sd ...
%                      + c_log_pair_dist + age + sex
%                      + (1 + delta_eeg_fc | subject_id) ...
%                      + (1 | subject_id:run_id) + (1 | pair_id)
%
% Notes:
% - Uses concatenated per-subject parts created by Step08A_build_windowLevel_LMM_parts...
% - Includes distance covariate (centered log distance).
% - Uses subject-specific random slope for delta_eeg_fc as the default model.
%
% RUN EXAMPLE:
% partsDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/parts_window";
% outDir   = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window";
% step08_fit_LMM_windowLevel_DEV_fromParts(partsDir, outDir);

% GROUP_MASTER/
% └── group_stats_outputs_window/
%     ├── LMM_window_delta_randomSlopeSubject_DEFAULT_fdrAcrossBands.csv
%     ├── LMM_window_delta_randomSlopeSubject_DEFAULT_fdrAcrossBands.mat
%     └── BLUPs/
%         ├── BLUP_deltaEEG_delta.csv
%         ├── BLUP_deltaEEG_theta.csv
%         ├── BLUP_deltaEEG_alpha.csv
%         ├── BLUP_deltaEEG_beta.csv
%         ├── BLUP_deltaEEG_low_gamma.csv
%         └── BLUP_deltaEEG_high_gamma.csv



if nargin < 1 || isempty(partsDir)
    partsDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/parts_window";
end
if nargin < 2 || isempty(outDir)
    outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window";
end
if ~exist(outDir,"dir"); mkdir(outDir); end

bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

% ---------- Load + concat all subjects ----------
files = dir(fullfile(partsDir, "ICE*_forLMM_window.mat"));
if isempty(files)
    error("No *_forLMM_window.mat files found in %s", partsDir);
end

Tall = table();
for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name), "Tsub");
    Tall = [Tall; S.Tsub]; %#ok<AGROW>
end

% ---------- Ensure required types ----------
Tall.subject_id = categorical(string(Tall.subject_id));

% % % ---------- OPTIONAL: remove outlier subjects ----------
% % Provide subject IDs as strings like "ICE013"
% removeSubj = ["ICE013","ICE028","ICE041","ICE050","ICE014"];  % <-- EDIT THIS LIST
% 
% removeSubj = categorical(removeSubj);
% nBefore = height(Tall);
% 
% isOut = ismember(Tall.subject_id, removeSubj);
% Tall(isOut,:) = [];
% 
% fprintf("[OUTLIER] Removed %d rows from %d subjects: %s\n", ...
%     sum(isOut), numel(removeSubj), strjoin(string(removeSubj), ", "));
% fprintf("[OUTLIER] Rows: %d -> %d (removed %.2f%%)\n", ...
%     nBefore, height(Tall), 100*(nBefore-height(Tall))/max(nBefore,1));
% 
% % ---------- END OF: remove outlier subjects ----------


Tall.pair_id    = categorical(string(Tall.pair_id));
Tall.band       = categorical(string(Tall.band), cellstr(bandOrder), "Ordinal", true);

if ismember("sex", Tall.Properties.VariableNames)
    Tall.sex = categorical(string(Tall.sex));
else
    Tall.sex = categorical(repmat("NA", height(Tall), 1));
end

% --- run-level random intercept support (REQUIRED) ---
if ~ismember("run_id", Tall.Properties.VariableNames)
    error("run_id not found — window-level model requires run information");
end
Tall.run_id = categorical(string(Tall.run_id));

% Drop rows with missing run_id
badRun = ismissing(Tall.run_id);
if any(badRun)
    warning("Found %d rows with missing run_id. Dropping them.", sum(badRun));
    Tall(badRun,:) = [];
end

% ---------- Distance covariate ----------
if ~ismember("pair_dist_mm", Tall.Properties.VariableNames)
    error("pair_dist_mm not found in Tall. Ensure Step07 + Step08A keep it and rebuild parts.");
end
Tall.pair_dist_mm = double(Tall.pair_dist_mm);

% Safe log transform (eps avoids log(0) just in case)
Tall.log_pair_dist = log(max(Tall.pair_dist_mm, eps));

% Center globally
Tall.c_log_pair_dist = Tall.log_pair_dist - mean(Tall.log_pair_dist, "omitnan");

% ---------- Single default model ----------
% formula = ...
%  "fmri_fc_windowed ~ eeg_fc_trueStatic + delta_eeg_fc + eeg_fc_windowed_sd + c_log_pair_dist + age + sex" + ...
%  " + (1|subject_id:run_id) + (1|pair_id) + (1 + delta_eeg_fc | subject_id)";

formula = "fmri_fc_windowed ~ eeg_fc_trueStatic + delta_eeg_fc + eeg_fc_windowed_sd" + ...
    " + c_log_pair_dist + delta_eeg_fc:c_log_pair_dist + age + sex" + ...
    " + (1|subject_id:run_id) + (1|pair_id)" + ...
    " + (1 + delta_eeg_fc + delta_eeg_fc:c_log_pair_dist | subject_id)";

% ---------- Fit per band ----------
rows = {};

blupOutDir = fullfile(outDir, "BLUPs");
if ~exist(blupOutDir,"dir"); mkdir(blupOutDir); end

for b = bandOrder
    Tb = Tall(Tall.band == b, :);

    useCols = ["fmri_fc_windowed","eeg_fc_trueStatic","delta_eeg_fc","eeg_fc_windowed_sd", ...
               "c_log_pair_dist","age","sex","subject_id","run_id","pair_id"];

    % Check missing columns
    missing = setdiff(string(useCols), string(Tb.Properties.VariableNames));
    nObs   = height(Tb);
    nSubj  = numel(unique(Tb.subject_id));
    nPairs = numel(unique(Tb.pair_id));


    if ~isempty(missing)
        fprintf("[SKIP] Band %s missing columns: %s\n", string(b), strjoin(missing, ", "));
        rows(end+1,:) = {b, ...
            NaN,NaN,NaN,NaN, ... % dev
            NaN,NaN,NaN,NaN, ... % static
            NaN,NaN,NaN,NaN, ... % sd
            NaN,NaN,NaN,NaN, ... % age
            NaN,NaN,NaN,NaN, ... % sex
            "", ...              % sex_coef_name
            false, false, "Missing columns", ...
            nObs_raw, nSubj_raw, nPairs_raw};
        continue;
    end

    % Drop rows with missing data in model variables
    Tb = rmmissing(Tb, "DataVariables", useCols);

    % Counts after cleaning
    nObs   = height(Tb);
    nSubj  = numel(unique(Tb.subject_id));
    nPairs = numel(unique(Tb.pair_id));


    if nObs < 1000 || nSubj < 2 || nPairs < 10
        rows(end+1,:) = {b, ...
            NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,NaN, ...
            "", ...
            false, false, "Too few observations/levels after cleaning", ...
            nObs, nSubj, nPairs};
        continue;
    end

    % --- Fit default model (keep model fit + row append separate from BLUP saving) ---
    converged = true;
    singular  = false;
    msg       = "";

    fprintf("\n[FIT] Band %s: nObs=%d  nSubj=%d  nPairs=%d\n", string(b), nObs, nSubj, nPairs);

    try
        lme = fitlme(Tb, formula, "FitMethod","REML");

        % singular check (version-dependent)
        try
            singular = lme.ModelCriterion.singular;
        catch
            singular = false;
        end

        % ---- Fixed effects (these should not depend on BLUP extraction) ----
        r_dev   = lmm_extract(lme, "delta_eeg_fc");
        r_stat  = lmm_extract(lme, "eeg_fc_trueStatic");
        r_sd    = lmm_extract(lme, "eeg_fc_windowed_sd");
        r_age   = lmm_extract(lme, "age");

        sexTerm = "";
        r_sex = struct("beta",NaN,"p",NaN,"ci_low",NaN,"ci_high",NaN);

        try
            sexTerm = find_sex_term(lme);
            r_sex   = lmm_extract(lme, sexTerm);
        catch MEsx
            % sex term may be dropped if only one level remains after filtering
            fprintf(2,"[WARN] Band %s: sex term not estimable (%s). Storing NaNs.\n", string(b), MEsx.message);
        end


        % ✅ Append ONE row only (never again in BLUP logic)
        rows(end+1,:) = {b, ...
            r_dev.beta,   r_dev.p,   r_dev.ci_low,   r_dev.ci_high, ...
            r_stat.beta,  r_stat.p,  r_stat.ci_low,  r_stat.ci_high, ...
            r_sd.beta,    r_sd.p,    r_sd.ci_low,    r_sd.ci_high, ...
            r_age.beta,   r_age.p,   r_age.ci_low,   r_age.ci_high, ...
            r_sex.beta,   r_sex.p,   r_sex.ci_low,   r_sex.ci_high, ...
            sexTerm, ...
            converged, singular, msg, ...
            nObs, nSubj, nPairs};

        % =================== BLUP extraction (safe / non-fatal) ===================
        try
            Tsubj = extract_subject_delta_blups(lme, "subject_id", "delta_eeg_fc", b, bandOrder);

            outCsv = fullfile(blupOutDir, sprintf("BLUP_deltaEEG_%s.csv", string(b)));
            writetable(Tsubj, outCsv);
            fprintf("[BLUP] Saved %s (nSubj=%d)\n", outCsv, height(Tsubj));

            % Print each subject slope (requested)
            for ii = 1:height(Tsubj)
                fprintf("[BLUP] %s %s (%s): slope=%.4f (dev=%.4f)\n", ...
                    string(b), string(Tsubj.subject_id(ii)), string(Tsubj.term(ii)), ...
                    Tsubj.deltaEEG_subjectSlope(ii), Tsubj.deltaEEG_subjectDeviation_from_group(ii));
            end

        catch MEblup
            fprintf(2, "[BLUP-WARN] Band %s: BLUP extraction failed: %s\n", string(b), MEblup.message);
        end

    catch ME
        converged = false;
        msg = ME.message;

        fprintf(2, "[FIT-FAIL] Band %s: %s\n", string(b), msg);

        % Append ONE failure row
        rows(end+1,:) = {b, ...
            NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,NaN, ...
            "", ...
            converged, singular, msg, ...
            nObs, nSubj, nPairs};
    end


end

% ---------- Convert to table ----------
Tout = cell2table(rows, "VariableNames", ...
    ["band", ...
    "beta_dev","p_dev","ci_low_dev","ci_high_dev", ...
    "beta_static","p_static","ci_low_static","ci_high_static", ...
    "beta_sd","p_sd","ci_low_sd","ci_high_sd", ...
    "beta_age","p_age","ci_low_age","ci_high_age", ...
     "beta_sex","p_sex","ci_low_sex","ci_high_sex", ...
     "sex_coef_name", ...
     "converged","singular","message", ...
     "nObs","nSubj","nPairs"]);

% ---------- BH-FDR across bands ----------
Tout.p_dev_fdr    = fdr_bh(Tout.p_dev);
Tout.p_static_fdr = fdr_bh(Tout.p_static);
Tout.p_sd_fdr     = fdr_bh(Tout.p_sd);
Tout.p_age_fdr    = fdr_bh(Tout.p_age);
Tout.p_sex_fdr    = fdr_bh(Tout.p_sex);

% ---------- Save outputs ----------
corrTag = "fdrAcrossBands";
base = "LMM_window_delta_randomSlopeSubject_DEFAULT_" + corrTag;

writetable(Tout, fullfile(outDir, base + ".csv"));
save(fullfile(outDir, base + ".mat"), "Tout", "-v7.3");

fprintf("[OK] Window-level DEV LMM complete (single DEFAULT model; BH-FDR across bands).\nSaved to: %s\n", outDir);

end

% ===================== Helpers =====================

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

function sexTerm = find_sex_term(lme)
names = string(lme.Coefficients.Name);
idx = startsWith(names, "sex_");

if sum(idx) == 1
    sexTerm = char(names(idx));
elseif sum(idx) == 0
    error("No sex coefficient found. Available terms: %s", strjoin(names, ", "));
else
    warning("Multiple sex coefficients found: %s. Using the first one.", strjoin(names(idx), ", "));
    sexTerm = char(names(find(idx,1,'first')));
end
end

function p_adj = fdr_bh(p)
p = double(p(:));
p(isnan(p)) = 1;

[ps, idx] = sort(p);
m = numel(p);
q = ps .* m ./ (1:m)';

q = flipud(cummin(flipud(q)));
q(q > 1) = 1;

p_adj = nan(size(p));
p_adj(idx) = q;
end

function Tsubj = extract_subject_delta_blups(lme, subjGroupName, slopeTerm, bandName, bandOrder)
% Robust extraction of subject-specific BLUP slopes across MATLAB versions.
% Tries new randomEffects(...,"Grouping",...) first; falls back to legacy outputs.

    subjLevels = strings(0,1);
    subjRand   = [];

    % ---------- Try NEW syntax (if supported) ----------
    try
        RE = randomEffects(lme, "Grouping", string(subjGroupName));  % newer MATLAB
        % Expect table with Level/Name/Estimate-like columns
        [subjLevels, subjRand, termName] = parse_RE_table(RE, slopeTerm);
    catch MEnew
        % If Grouping isn't supported, fall back
        if ~contains(MEnew.message, "Grouping", "IgnoreCase", true)
            % Grouping failed for another reason; still fall back anyway
        end

        % ---------- Try legacy table output ----------
        try
            RE1 = randomEffects(lme);  % older MATLAB: may return table OR numeric
            if istable(RE1) || isa(RE1,"dataset") || isstruct(RE1)
                if isa(RE1,"dataset"), RE1 = dataset2table(RE1); end
                if isstruct(RE1),      RE1 = struct2table(RE1);  end

                % Legacy tables usually have Group/Level/Name/Estimate
                [subjLevels, subjRand, termName] = parse_RE_table_legacy(RE1, subjGroupName, slopeTerm);
            else
                error("randomEffects(lme) not table-like; switching to [b,names] mode.");
            end

        catch
            % ---------- Last resort: numeric + names ----------
            [b, names] = randomEffects(lme);
            b = double(b(:));

            % ---- CASE: names is a TABLE (common in older MATLAB) ----
            if istable(names)
                vn = string(names.Properties.VariableNames);
                low = lower(vn);

                % Try to find the key columns
                groupVar = vn(find(contains(low,"group"), 1, "first"));
                levelVar = vn(find(contains(low,"level"), 1, "first"));
                nameVar  = vn(find(contains(low,"name") | contains(low,"term"), 1, "first"));

                if isempty(groupVar) || isempty(levelVar) || isempty(nameVar)
                    error("randomEffects second output is table but missing Group/Level/Name columns. Have: %s", strjoin(vn,", "));
                end

                groupCol = string(names.(groupVar));
                levelCol = string(names.(levelVar));
                nameCol  = string(names.(nameVar));

                % keep subject group only
                isSubj = (groupCol == string(subjGroupName)) | contains(groupCol, string(subjGroupName));
                if ~any(isSubj)
                    error("No subject-level random effects found in names table. Groups: %s", strjoin(unique(groupCol), ", "));
                end

                % keep slopeTerm only
                levelCol = levelCol(isSubj);
                nameCol  = nameCol(isSubj);
                b2       = b(isSubj);

                isSlope = (nameCol == string(slopeTerm)) | startsWith(nameCol, string(slopeTerm));
                if ~any(isSlope)
                    error("No random slope '%s' found. Names: %s", slopeTerm, strjoin(unique(nameCol), ", "));
                end

                termName = nameCol(isSlope);
                subjLevels = levelCol(isSlope);
                subjRand   = b2(isSlope);

            else
                % ---- CASE: names is cell/string/char vector ----
                if iscell(names), names = string(names); end
                names = string(names(:));

                % DEBUG safely
                fprintf("[DEBUG] randomEffects names example:\n");
                disp(names(1:min(10,height(names)), :));

                [subjLevels, subjRand, termName] = parse_RE_names(names, b, subjGroupName, slopeTerm);

            end

        end
    end

    if isempty(subjLevels)
        error("Could not extract subject random slopes for '%s' from randomEffects outputs.", slopeTerm);
    end

    % ---- Fixed slope + random deviation = subject-specific slope ----
    fixedSlope = lmm_extract(lme, slopeTerm).beta;

    Tsubj = table();
    Tsubj.subject_id                           = categorical(subjLevels);
    Tsubj.deltaEEG_subjectDeviation_from_group = subjRand(:);
    Tsubj.deltaEEG_subjectSlope                = fixedSlope + subjRand(:);
    Tsubj.band                                 = categorical(repmat(string(bandName), numel(subjRand), 1), cellstr(bandOrder), "Ordinal", true);
    Tsubj.term                                 = categorical(termName);

end

% ===== Helpers for parsing randomEffects outputs =====

function [levels, est, termName] = parse_RE_table(RE, slopeTerm)
    v = string(RE.Properties.VariableNames);
    nameVar = v(contains(lower(v),"name", "IgnoreCase", true));  nameVar = nameVar(1);
    levVar  = v(contains(lower(v),"level","IgnoreCase", true));  levVar  = levVar(1);
    estVar  = v(contains(lower(v),"estimate","IgnoreCase", true) | contains(lower(v),"est","IgnoreCase", true)); estVar = estVar(1);

    name = string(RE.(nameVar));
    lev  = string(RE.(levVar));
    est0 = double(RE.(estVar));

    isSlope = (name == string(slopeTerm)) | startsWith(name, string(slopeTerm));

    levels   = lev(isSlope);
    est      = est0(isSlope);
    termName = name(isSlope);   % <<< NEW
end


function [levels, est, termName] = parse_RE_table_legacy(RE, subjGroupName, slopeTerm)
    v = string(RE.Properties.VariableNames);
    low = lower(v);

    groupVar = v(find(contains(low,"group") | contains(low,"grouping"), 1, "first"));
    levVar   = v(find(contains(low,"level"), 1, "first"));
    nameVar  = v(find(contains(low,"name") | contains(low,"term"), 1, "first"));
    estVar   = v(find(contains(low,"estimate") | contains(low,"est"), 1, "first"));

    groupCol = string(RE.(groupVar));
    levCol   = string(RE.(levVar));
    nameCol  = string(RE.(nameVar));
    estCol   = double(RE.(estVar));

    isSubj = (groupCol == string(subjGroupName)) | contains(groupCol, string(subjGroupName));
    if ~any(isSubj)
        error("No subject-level random effects found. Group values: %s", strjoin(unique(groupCol), ", "));
    end

    levCol  = levCol(isSubj);
    nameCol = nameCol(isSubj);
    estCol  = estCol(isSubj);

    isSlope = (nameCol == string(slopeTerm)) | startsWith(nameCol, string(slopeTerm));
    if ~any(isSlope)
        error("No random slope '%s' found. Names: %s", slopeTerm, strjoin(unique(nameCol), ", "));
    end

    levels   = levCol(isSlope);
    est      = estCol(isSlope);
    termName = nameCol(isSlope);   % <<< NEW
end


function [levels, est, termName] = parse_RE_names(names, b, subjGroupName, slopeTerm)
    isSlope = contains(names, "):" + slopeTerm) & contains(names, subjGroupName + "(");

    candNames = names(isSlope);
    candB     = b(isSlope);

    if isempty(candNames)
        error("Could not find any slope entries for %s in randomEffects names.", slopeTerm);
    end

    levels = strings(numel(candNames),1);
    for i = 1:numel(candNames)
        tmp = extractBetween(candNames(i), subjGroupName+"(", "):");
        if isempty(tmp)
            tmp = extractBetween(candNames(i), "(", ")");
        end
        levels(i) = string(tmp);
    end

    est      = candB(:);
    termName = candNames;  % <<< NEW: store full name string so you know which term it is
end



function vname = pick_var(vars, low, candidates)
% pick the first matching column name among candidates
    vname = "";
    for c = candidates
        idx = (low == lower(string(c)));
        if any(idx)
            vname = vars(find(idx,1,'first'));
            return;
        end
    end
    % try fuzzy contains
    for c = candidates
        idx = contains(low, lower(string(c)));
        if any(idx)
            vname = vars(find(idx,1,'first'));
            return;
        end
    end
    error("Could not find required column. Have columns: %s", strjoin(vars, ", "));
end

