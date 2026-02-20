function step08A_build_windowLevel_LMM_parts_fromMasterTables_v2()
% prepares the window-level analysis table (“LMM-ready”) by joining:
% per-window rows (from master_dynamicFC.csv)
% with
% pair-level predictors (from master_staticFC.csv)
% and computing within-subject deviation:
% delta_eeg_fc is already computed in Step07:
% delta_eeg_fc = eeg_fc_windowed - eeg_fc_trueStatic
% delta_eeg_fc = within-subject, within-pair EEG FC deviation
%                (windowed EEG FC minus subject-level true static EEG FC)

% Builds per-subject, LMM-ready WINDOW tables and saves as .mat parts.

% IMPORTANT: keeps pair_id, run_id, win_idx so we can fit:
%   (1 | subject_id) + (1 | subject_id:run_id) + (1 | pair_id)

%
% Output:
%   .../master_tables/GROUP_MASTER/parts/ICE###_forLMM_window.mat

baseDir  = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables";
subjDirs = dir(fullfile(baseDir,"ICE*"));
subjDirs = subjDirs([subjDirs.isdir]);

groupDir = fullfile(baseDir,"GROUP_MASTER");
partsDir = fullfile(groupDir,"parts_window");
if ~exist(partsDir,"dir"); mkdir(partsDir); end

bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

for i = 1:numel(subjDirs)
    subj = string(subjDirs(i).name);

    winPath  = fullfile(baseDir, subj, "master_dynamicFC.csv");
    pairPath = fullfile(baseDir, subj, "master_staticFC.csv");

    if ~exist(winPath,"file") || ~exist(pairPath,"file")
        fprintf("[SKIP] %s missing master tables\n", subj);
        continue;
    end

    Tw = readtable(winPath);
    Tp = readtable(pairPath);

    % Keep only the columns we need from pair table
    Tp2 = Tp(:, {'subject_id','pair_id','band','eeg_fc_trueStatic','eeg_fc_windowed_sd'});

    % Normalize join key types
    Tw.subject_id = string(Tw.subject_id);
    Tw.pair_id    = string(Tw.pair_id);
    Tw.band       = string(Tw.band);
    Tw.run_id     = string(Tw.run_id);

    Tp2.subject_id = string(Tp2.subject_id);
    Tp2.pair_id    = string(Tp2.pair_id);
    Tp2.band       = string(Tp2.band);

    % Join to bring eeg_fc_static + eeg_fc_dynamic_sd into window rows
    Tj = innerjoin(Tw, Tp2, 'Keys', {'subject_id','pair_id','band'});

    % ---- SAFETY: ensure pair_dist_mm exists after join ----
    if ~ismember("pair_dist_mm", Tj.Properties.VariableNames)
        % if join created suffixed versions, pick one
        cand = Tj.Properties.VariableNames(contains(Tj.Properties.VariableNames, "pair_dist_mm"));
        error("pair_dist_mm missing after join. Found instead: %s", strjoin(string(cand), ", "));
    end

    % Ensure numeric
    Tj.pair_dist_mm = double(Tj.pair_dist_mm);

    % Recommended transforms:
    % Robust log (avoid log(0))
    Tj.log_pair_dist = log(Tj.pair_dist_mm + eps);           % distance is often right-skewed
    Tj.c_log_pair_dist = Tj.log_pair_dist - mean(Tj.log_pair_dist, 'omitnan');  % center


    if height(Tj) ~= height(Tw)
        fprintf("[WARN] %s join dropped rows: Tw=%d -> Tj=%d (%.1f%% retained)\n", ...
            subj, height(Tw), height(Tj), 100*height(Tj)/height(Tw));
    end

    % --- compute within-subject deviation term (old manuscript) ---
%     Tj.delta_eeg_fc = Tj.eeg_fc_windowed - Tj.eeg_fc_trueStatic;

    % Standardize band naming/order
    Tj.band = string(Tj.band);
    Tj.band(Tj.band=="gammalow")  = "low_gamma";
    Tj.band(Tj.band=="gammahigh") = "high_gamma";
    Tj.band = categorical(Tj.band, cellstr(bandOrder), "Ordinal", true);

    % standardize subject_id and sex for fitlme
    Tj.subject_id = categorical(string(Tj.subject_id));
    Tj.pair_id    = categorical(string(Tj.pair_id));
    Tj.run_id     = categorical(string(Tj.run_id));
    
    % Ensure age/sex exist (some subjects may miss these)
    if ismember("sex", Tj.Properties.VariableNames)
        Tj.sex = categorical(string(Tj.sex));
    else
        Tj.sex = categorical(repmat("NA", height(Tj), 1));
    end

    if ~ismember("age", Tj.Properties.VariableNames)
        Tj.age = nan(height(Tj),1);    
    end

    % Ensure pair-level predictors exist (should after join, but guard anyway)
    if ~ismember("eeg_fc_trueStatic", Tj.Properties.VariableNames)
        error("%s is missing eeg_fc_trueStatic after join. Check master_staticFC.csv", subj);
    end
    if ~ismember("eeg_fc_windowed_sd", Tj.Properties.VariableNames)
        error("%s is missing eeg_fc_windowed_sd after join. Check master_staticFC.csv", subj);
    end

    if ~ismember("delta_eeg_fc", Tj.Properties.VariableNames)
        error("%s is missing delta_eeg_fc. Check master_dynamicFC.csv creation in Step07.", subj);
    end
    if ~ismember("fmri_fc_windowed", Tj.Properties.VariableNames)
        error("%s is missing fmri_fc_windowed. Check master_dynamicFC.csv.", subj);
    end


    % Keep essential columns (KEEP pair/run/win)

    keep = {'subject_id','pair_id','band','run_id','win_idx', ...
    'fmri_fc_windowed','eeg_fc_windowed','delta_eeg_fc', ...
    'eeg_fc_trueStatic','eeg_fc_windowed_sd', ...
    'pair_dist_mm','log_pair_dist','c_log_pair_dist', ...
    'age','sex'};


    missing = setdiff(keep, Tj.Properties.VariableNames);
    if ~isempty(missing)
        error("%s missing required window-level columns: %s", subj, strjoin(missing, ", "));
    end

    % Some older files may not have all columns; guard:
    keep = keep(ismember(keep, Tj.Properties.VariableNames));
    Tsub = Tj(:, keep);

    outMat = fullfile(partsDir, subj + "_forLMM_window.mat");
    Tsub = sortrows(Tsub, {'subject_id','pair_id','run_id','win_idx'});
    save(outMat, "Tsub", "-v7.3");

    fprintf("[OK] Saved %s rows=%d -> %s\n", subj, height(Tsub), outMat);
end

fprintf("\n[DONE] Window-level LMM parts saved to:\n%s\n", partsDir);
end
