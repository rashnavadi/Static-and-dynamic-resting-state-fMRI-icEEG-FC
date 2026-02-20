function Step07_build_master_tables(cfgPaths, masterOutDir, demoCsvPath, makeWindowTable)
% master_dynamicFC.csv: *_dynamic = per-window
% master_staticFC.csv: *_dynamic_mean = mean across windows/runs

% Step07 computes:
% midpoint of each bipolar channel (contact A/B midpoint)
% Euclidean distance between the two bipolar midpoints for each edge/pair
% stored as pair_dist_mm (length = nPairs)

% Step07_build_master_tables
% Builds frozen "master" tables across subjects from:
% Inputs:
%   - Step04: FMRI_FC_<subj>.mat  (OUT struct)
%   - Step05: EEG_FC_Hilbert_<subj>.mat (OUT struct)
%
% Outputs:
%   master_staticFC.csv/.mat   (subject × pair × band)
%   master_dynamicFC.csv/.mat (subject × pair × band × window) [optional]
%   bipolar midpoint CSV we just added goes: outBipCsv = fullfile(coordsDir, subj + "_bipolar_nodes_midpoints.csv");

%
% IMPORTANT:
% - Uses Fisher-z consistently for BOTH modalities.
% - Maps gammalow->low_gamma, gammahigh->high_gamma.
% - Dynamic (pair-level) is MEAN across ALL windows and ALL runs per subject.
% - SD is computed across ALL windows and ALL runs per subject (on Fisher-z values).

if nargin < 4
    makeWindowTable = true; % set false if you want only pair-level
end

if ~exist(masterOutDir, "dir")
    mkdir(masterOutDir);
end

% ---------- Load electrodes coords in native space ----------

coordsDir = "/Volumes/MIND/ICE/Tara/native_space_electrodes_coords";

% ---------- Load demographics/covariates ----------
demo = readtable(demoCsvPath);
% force exact column name
if ~ismember("subject_id", string(demo.Properties.VariableNames))
    error("demographics file must contain column named 'subject_id'");
end
if ~ismember("age", string(demo.Properties.VariableNames))
    error("demographics file must contain column named 'age'");
end
if ~ismember("sex", string(demo.Properties.VariableNames))
    error("demographics file must contain column named 'sex'");
end

demo.subject_id = string(demo.subject_id);
demo.sex        = string(demo.sex);

% Manuscript band order
bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

subjSummaryRows = {};  % will become a table at end

for s = 1:numel(cfgPaths)
    cfgPath = string(cfgPaths{s});
    S = load(cfgPath, "C");
    C = S.C;

    % Ensure helpers are on path
    if isfield(C, "codeDir")
        addpath(C.codeDir);
        addpath(fullfile(C.codeDir, "helpers"));
    end

    subj = string(C.subjectID);
    fprintf("\n=== Building master rows for %s ===\n", subj);

    clear pair_dist_mm

    % create a subject subfolder
    subjOutDir = fullfile(masterOutDir, subj);
    if ~exist(subjOutDir, "dir")
        mkdir(subjOutDir);
    end

    % Load Step04 and Step05 outputs
    fmriFCpath = fullfile(C.outDir, sprintf("FMRI_FC_%s.mat", subj));
    eegFCpath  = fullfile(C.outDir, sprintf("EEG_FC_Hilbert_%s.mat", subj));

    if ~exist(fmriFCpath,"file"); error("Missing %s", fmriFCpath); end
    if ~exist(eegFCpath,"file");  error("Missing %s", eegFCpath);  end

    F = load(fmriFCpath, "OUT");
    E = load(eegFCpath,  "OUT");

    OUTf = F.OUT;
    OUTe = E.OUT;

    % ---------- Covariates from demographics.csv ----------
    idx = find(demo.subject_id == subj, 1);
    if isempty(idx)
        error("Subject %s not found in demographics CSV: %s", subj, demoCsvPath);
    end
    age = demo.age(idx);
    sex = demo.sex(idx);   % already string

    % Accumulators across runs (subject-level)
    subj_run_static_fmri = [];  % [nPairs x nRuns]
    subj_run_fmri_allWindowed = []; % [nPairs x totalWins]
    % EEG accumulators per band
    eeg_static_byBand = struct();
    eeg_dyn_allWins_byBand = struct();
    for b = bandOrder
        bn = char(b);
        eeg_static_byBand.(bn) = [];
        eeg_dyn_allWins_byBand.(bn) = [];
    end

    % For window table global win indexing within subject
    globalWinOffset = 0;

    pairRows_subj = {};
    winRows_subj  = {};

    runKeys = intersect(fieldnames(OUTf), fieldnames(OUTe), "stable");
    if isempty(runKeys)
        warning("No common run keys between fMRI and EEG OUT structs for %s. Skipping subject.", subj);
        continue;
    end
    pairNames_subject = strings(0,1);
    nPairs_subject = [];
    roiCounts_validRuns = [];   % store N for each valid run used
    nValidRuns = 0;

    for rk = 1:numel(runKeys)
        runKey = runKeys{rk};
        if ~isstruct(OUTf.(runKey)) || ~isstruct(OUTe.(runKey))
            continue;
        end

        % ---------- Labels / names ----------
        names_fmri = string(OUTf.(runKey).seedNames);
        names_fmri = upper(strrep(names_fmri,"_","-"));
        names_fmri = canonical_bipolar_name(names_fmri);

        names_eeg  = string(OUTe.(runKey).labels);
        names_eeg  = upper(strrep(names_eeg,"_","-"));
        names_eeg  = canonical_bipolar_name(names_eeg);

        % If EEG also has pairKeepVec, enforce consistency
        if isfield(OUTe.(runKey), "pairKeepVec") && isfield(OUTf.(runKey), "pairKeepVec")
            keepVec_e = logical(OUTe.(runKey).pairKeepVec(:));
            keepVec_f = logical(OUTf.(runKey).pairKeepVec(:));
            if numel(keepVec_e) == numel(keepVec_f) && any(keepVec_e ~= keepVec_f)
                warning("EEG/FMRI pairKeepVec mismatch in %s %s. Skipping run.", subj, runKey);
                continue;
            end
        end

        % Must match order exactly (your Step05 enforces FMRI order)
        if numel(names_fmri) ~= numel(names_eeg) || any(names_fmri ~= names_eeg)
            warning("Name/order mismatch in %s %s. Skipping this run.", subj, runKey);
            continue;
        end

        N = numel(names_fmri);

        % Full upper-triangle pair names (length = N*(N-1)/2)
        pairNames_full = uppertri_pairnames(names_fmri);

        if ~isfield(OUTf.(runKey), "pairKeepVec") || isempty(OUTf.(runKey).pairKeepVec)
            warning("Missing pairKeepVec in %s %s. Skipping run to avoid mixing full vs kept edges.", subj, runKey);
            continue;
        end

        % Use survived-edge mask if present (preferred)
        if isfield(OUTf.(runKey), "pairKeepVec") && ~isempty(OUTf.(runKey).pairKeepVec)
            keepVec_f = logical(OUTf.(runKey).pairKeepVec(:));

            if numel(keepVec_f) ~= numel(pairNames_full)
                warning("pairKeepVec length mismatch in %s %s. Expected %d got %d. Skipping run.", ...
                    subj, runKey, numel(pairNames_full), numel(keepVec_f));
                continue;
            end

            pairNames = pairNames_full(keepVec_f);
            nPairs = numel(pairNames);
        else
            % Fallback to full graph
            pairNames = pairNames_full;
            nPairs = numel(pairNames_full);
        end

        % Freeze pairNames at subject level on first valid run
        if isempty(nPairs_subject)
            nPairs_subject = nPairs;
            pairNames_subject = pairNames;
        else
            % Enforce same pair ordering/length across runs used for accumulation
            if nPairs ~= nPairs_subject || any(pairNames ~= pairNames_subject)
                warning("Pair list mismatch across runs in %s. Skipping run %s.", subj, runKey);
                continue;
            end
        end

        fprintf("[PAIRSET] %s %s: N=%d fullPairs=%d keptPairs=%d\n", ...
            subj, runKey, N, numel(pairNames_full), nPairs);

        if isfield(OUTf.(runKey), "nPairsKeep")
            fprintf("          OUTf nPairsKeep=%d\n", OUTf.(runKey).nPairsKeep);
        end
        if isfield(OUTe.(runKey), "nPairsKeep")
            fprintf("          OUTe nPairsKeep=%d\n", OUTe.(runKey).nPairsKeep);
        end


        % ================== DISTANCE COMPUTATION (ONCE PER SUBJECT) ==================
        if ~exist("pair_dist_mm", "var") || isempty(pair_dist_mm)
            pair_dist_mm = compute_pair_distance_mm(subj, pairNames_subject, coordsDir);
        end

        % SAFETY CHECK (place it HERE, right after distance is created)
        if numel(pair_dist_mm) ~= numel(pairNames_subject)
            error("[DIST] %s: pair_dist_mm length (%d) != pairNames_subject (%d)", ...
                subj, numel(pair_dist_mm), numel(pairNames_subject));
        end

        % ============================================================================
        % ---------- DEBUG + SAFETY: FC_static_z must match matched seed list ----------
        % ---------- fMRI static (Fisher-z matrix exists) ----------
        % Prefer survived-edge static vector; only require FC_static_z if staticVec_z is missing.
        if isfield(OUTf.(runKey), "staticVec_z") && ~isempty(OUTf.(runKey).staticVec_z)
            fmri_static_vec_z = OUTf.(runKey).staticVec_z(:);
        else
            if ~isfield(OUTf.(runKey), "FC_static_z") || isempty(OUTf.(runKey).FC_static_z)
                warning("Missing both staticVec_z and FC_static_z in %s %s. Skipping run.", subj, runKey);
                continue;
            end
            fmri_static_vec_z = uppertri_vec(OUTf.(runKey).FC_static_z);
        end


        % Prefer survived-edge static vector if present
        if isfield(OUTf.(runKey), "staticVec_z") && ~isempty(OUTf.(runKey).staticVec_z)
            fmri_static_vec_z = OUTf.(runKey).staticVec_z(:);
        else
            fmri_static_vec_z = uppertri_vec(OUTf.(runKey).FC_static_z);
        end

        if numel(fmri_static_vec_z) ~= nPairs
            fprintf("[DEBUG] %s %s: expected nPairs=%d but numel(fmri_static_vec_z)=%d\n", ...
                subj, runKey, nPairs, numel(fmri_static_vec_z));
            warning("Static fMRI vector length mismatch in %s %s. Skipping run.", subj, runKey);
            continue;
        end

        % Initialize subj_run_static_fmri if it's the first valid run
        if isempty(subj_run_static_fmri)
            subj_run_static_fmri = nan(nPairs, 0);
        end

        % ---------- fMRI dynamic (use dynVec_z directly) ----------
        fmri_dyn_z = OUTf.(runKey).dynVec_z; % [nPairs x nWin_fmri]
        if isempty(fmri_dyn_z)
            warning("Empty fMRI dynVec_z in %s %s. Skipping run.", subj, runKey);
            continue;
        end
        if size(fmri_dyn_z,1) ~= nPairs
            fprintf("[DEBUG] %s %s: nPairs(from names)=%d but size(dynVec_z,1)=%d\n", ...
                subj, runKey, nPairs, size(fmri_dyn_z,1));
            warning("fMRI dynVec_z row mismatch in %s %s. Skipping run.", subj, runKey);
            continue;
        end

        nWin_fmri = size(fmri_dyn_z, 2);

        % ---------- EEG: per band ----------
        % Map Step05 band fieldnames to manuscript band names
        bandMap = struct( ...
            "delta","delta", ...
            "theta","theta", ...
            "alpha","alpha", ...
            "beta","beta", ...
            "gammalow","low_gamma", ...
            "gammahigh","high_gamma" ...
            );

        % Find EEG bands available in this run
        eegBandFields = fieldnames(OUTe.(runKey));
        % Keep only those that are in bandMap
        eegBandFields = eegBandFields(ismember(eegBandFields, fieldnames(bandMap)));

        % Determine nWin_eeg from first available band
        nWin_eeg = [];
        if ~isempty(eegBandFields)
            firstBn = eegBandFields{1};
            tmpDyn = OUTe.(runKey).(firstBn).dynVec_z;
            if ~isempty(tmpDyn)
                nWin_eeg = size(tmpDyn,2);
            end
        end
        if isempty(nWin_eeg)
            warning("Cannot determine EEG window count in %s %s. Skipping run.", subj, runKey);
            continue;
        end

        % Align windows (truncate to common)
        nWin = min(nWin_fmri, nWin_eeg);
        if nWin <= 0
            warning("No overlapping windows in %s %s.", subj, runKey);
            continue;
        end
        subj_run_static_fmri(:, end+1) = fmri_static_vec_z;
        
        roiCounts_validRuns(end+1) = N;
        nValidRuns = nValidRuns + 1;

        % Append fMRI dynamic windows (all bands share same fMRI)
        subj_run_fmri_allWindowed = [subj_run_fmri_allWindowed, fmri_dyn_z(:, 1:nWin)];

        % Build window-level rows (optional) — do it per band to match manuscript figures
        % We need EEG dyn for each band; fMRI dyn is shared.
        for b = 1:numel(eegBandFields)
            bn_step05 = eegBandFields{b};                 % e.g., gammalow
            bn_master = char(bandMap.(bn_step05));      % low_gamma

            % Fisher-z EEG
            eeg_static_z = OUTe.(runKey).(bn_step05).staticVec_z;   % already z
            if numel(eeg_static_z) ~= nPairs
                fprintf("[DEBUG] %s %s %s: EEG nPairs=%d but numel(eeg_static_z)=%d\n", ...
                    subj, runKey, bn_master, nPairs, numel(eeg_static_z));
                warning("EEG staticVec_z length mismatch in %s %s (%s). Skipping band.", subj, runKey, bn_master);
                continue;
            end

            eeg_dyn_z    = OUTe.(runKey).(bn_step05).dynVec_z(:,1:nWin);
            if size(eeg_dyn_z,1) ~= nPairs
                fprintf("[DEBUG] %s %s %s: EEG nPairs=%d but eeg_dyn_z rows=%d\n", ...
                    subj, runKey, bn_master, nPairs, size(eeg_dyn_z,1));
                warning("EEG dynVec_z row mismatch in %s %s (%s). Skipping band.", subj, runKey, bn_master);
                continue;
            end

            % Compute delta EEG (within-subject deviation)
            delta_eeg_z = eeg_dyn_z - eeg_static_z;


            % Accumulate subject-level across runs
            eeg_static_byBand.(bn_master) = [eeg_static_byBand.(bn_master), eeg_static_z];
            eeg_dyn_allWins_byBand.(bn_master) = [eeg_dyn_allWins_byBand.(bn_master), eeg_dyn_z];

            % Optional window-level table rows for figures
            if makeWindowTable
                % Global window index within subject
                win_idx = (1:nWin) + globalWinOffset;

                % Create rows: each pair × each window
                % (store as cells for fast append, convert to table at end)
                fmri_dyn_z_trunc = fmri_dyn_z(:,1:nWin);

                for w = 1:nWin
                    % For this window w: vectors over pairs
                    fmri_w = fmri_dyn_z_trunc(:, w);     % already Fisher-z
                    eeg_w  = eeg_dyn_z(:, w);      % Fisher-z
                    delta_w = delta_eeg_z(:, w);

                    % Append nPairs rows
                    winRows_subj = [winRows_subj; ...
                        [repmat({subj}, nPairs,1), ...
                        repmat({string(runKey)}, nPairs,1), ...
                        cellstr(pairNames_subject), ...
                        num2cell(pair_dist_mm), ...
                        repmat({char(bn_master)}, nPairs,1), ...
                        num2cell(repmat(win_idx(w), nPairs,1)), ...
                        num2cell(fmri_w), num2cell(eeg_w), num2cell(delta_w), ...
                        num2cell(repmat(age, nPairs,1)), repmat({char(sex)}, nPairs,1)] ...
                        ];

                end
            end
        end

        globalWinOffset = globalWinOffset + nWin;
        fprintf("[OK] %s %s: N=%d pairs=%d windows=%d (aligned)\n", subj, runKey, N, nPairs, nWin);
    end


    % ---------- Subject-level finalize: fMRI ----------
    % TERMINOLOGY:
    % trueStatic = computed on full run (Step04 FC_static_z / Step05 staticVec_z)
    % windowed   = per-window FC (Step04 dynVec_z / Step05 dynVec_z)
    % windowAvg  = mean(windowed) across all windows & runs  (static summary used for comparability)

    if isempty(subj_run_static_fmri) || isempty(subj_run_fmri_allWindowed)
        warning("No valid runs accumulated for subject %s. Skipping subject.", subj);
        continue;
    end

    % Average static across runs (Fisher-z)
    fmri_static_mean_z = mean(subj_run_static_fmri, 2, "omitnan"); % [nPairs x 1]
    % Window-averaged fMRI FC (mean over windowed FC across all runs/windows)
    fmri_windowAvg_z = mean(subj_run_fmri_allWindowed, 2, "omitnan"); % [nPairs x 1]

    % We need pairNames; reuse from last valid run (safe because ordering fixed by names)
    % If you want, you can store pairNames once when first valid run found.

    % ---------- Subject-level finalize: EEG per band ----------
    for b = bandOrder
        bn = char(b);

        if isempty(eeg_static_byBand.(bn)) || isempty(eeg_dyn_allWins_byBand.(bn))
            % band not available for this subject
            continue;
        end

        eeg_static_mean_z = mean(eeg_static_byBand.(bn), 2, "omitnan");        % [nPairs x 1]
        % Window-averaged EEG FC (mean over windowed FC across all runs/windows)
        eeg_windowAvg_z = mean(eeg_dyn_allWins_byBand.(bn), 2, "omitnan"); % [nPairs x 1]
        eeg_dyn_sd_z      = std(eeg_dyn_allWins_byBand.(bn), 0, 2, "omitnan"); % [nPairs x 1] across all windows

        % Append pair-level rows: subject × pair × band
        nPairs = numel(eeg_static_mean_z);

        % If pairNames was not defined (edge case), derive from last run’s names
        if isempty(pairNames_subject)
            warning("pairNames_subject missing for %s; cannot write pair-level rows.", subj);
            continue;
        end

        pairRows_subj = [pairRows_subj; ...
            [repmat({subj}, nPairs,1), ...
            cellstr(pairNames_subject), ...
            num2cell(pair_dist_mm), ...
            repmat({bn}, nPairs,1), ...
            num2cell(fmri_static_mean_z), num2cell(fmri_windowAvg_z), ...
            num2cell(eeg_static_mean_z),  num2cell(eeg_windowAvg_z), num2cell(eeg_dyn_sd_z), ...
            num2cell(repmat(age, nPairs,1)), repmat({char(sex)}, nPairs,1)] ];

    end

    fprintf("[DONE] Subject %s appended to master rows.\n", subj);

    if isempty(roiCounts_validRuns)
        warning("No valid ROI counts for %s; skipping summary row.", subj);
    else
        roi_min = min(roiCounts_validRuns);
        roi_max = max(roiCounts_validRuns);

        % Choose one subject-level ROI number for reporting:
        % Option A (common): use roi_max as "final included ROIs"
        roi_subject = roi_max;

        % nPairs_subject is already tracked when runs are valid
        nPairs_subject_final = nPairs_subject;

        % total aligned windows used across all valid runs for this subject:
        nWins_total = globalWinOffset;  % because you add nWin each run

        subjSummaryRows = [subjSummaryRows; ...
            {char(subj), age, char(sex), roi_min, roi_max, roi_subject, ...
            nPairs_subject_final, nValidRuns, nWins_total}];
    end


    % ---------- Convert to tables ----------
    % ---------- Write per-subject outputs ----------
    Tpair_subj = cell2table(pairRows_subj, "VariableNames", ...
        ["subject_id","pair_id","pair_dist_mm","band", ...
        "fmri_fc_trueStatic","fmri_fc_windowAvg", ...
        "eeg_fc_trueStatic","eeg_fc_windowAvg","eeg_fc_windowed_sd", ...
        "age","sex"]);

    Tpair_subj.band = categorical(string(Tpair_subj.band), cellstr(bandOrder), "Ordinal", true);
    Tpair_subj.subject_id = categorical(string(Tpair_subj.subject_id));
    Tpair_subj.pair_id    = categorical(string(Tpair_subj.pair_id));
    Tpair_subj.sex        = categorical(string(Tpair_subj.sex));

    writetable(Tpair_subj, fullfile(subjOutDir, "master_staticFC.csv"));
    save(fullfile(subjOutDir, "master_staticFC.mat"), "Tpair_subj", "-v7.3");

    if makeWindowTable && ~isempty(winRows_subj)
        Twin_subj = cell2table(winRows_subj, "VariableNames", ...
            ["subject_id","run_id","pair_id","pair_dist_mm","band","win_idx", ...
            "fmri_fc_windowed","eeg_fc_windowed","delta_eeg_fc","age","sex"]);


        Twin_subj.band = categorical(string(Twin_subj.band), cellstr(bandOrder), "Ordinal", true);
        Twin_subj.subject_id = categorical(string(Twin_subj.subject_id));
        Twin_subj.run_id     = categorical(string(Twin_subj.run_id));
        Twin_subj.pair_id    = categorical(string(Twin_subj.pair_id));
        Twin_subj.sex        = categorical(string(Twin_subj.sex));

        writetable(Twin_subj, fullfile(subjOutDir, "master_dynamicFC.csv"));
        save(fullfile(subjOutDir, "master_dynamicFC.mat"), "Twin_subj", "-v7.3");

    end
    fprintf("[OK] Saved per-subject master tables to %s\n", subjOutDir);

end
if ~isempty(subjSummaryRows)
    TsubjSummary = cell2table(subjSummaryRows, "VariableNames", ...
        ["subject_id","age","sex","roi_min","roi_max","roi_subject", ...
        "n_pairs","n_valid_runs","n_windows_total"]);

    writetable(TsubjSummary, fullfile(masterOutDir, "master_subject_roi_pair_summary.csv"));
    save(fullfile(masterOutDir, "master_subject_roi_pair_summary.mat"), "TsubjSummary");
end


end

% -------------------- Helpers (local) --------------------

function pairNames = uppertri_pairnames(names)
% names: [N x 1] string
names = string(names(:));
N = numel(names);
pairNames = strings(N*(N-1)/2, 1);
k = 1;
for i = 1:N-1
    for j = i+1:N
        pairNames(k) = names(i) + "--" + names(j);
        k = k + 1;
    end
end
end

function mid = get_midpoint_for_bipolar(chanName, contactMap, chanMidCache)
% chanName: e.g., "DRA1-DRA2" (already canonicalized upstream)
% contactMap: containers.Map contactKey-> [x y z]
% chanMidCache: containers.Map for caching midpoints by chanName

key = char(chanName);
if isKey(chanMidCache, key)
    mid = chanMidCache(key);
    return;
end

parts = strsplit(key, "-");
if numel(parts) ~= 2
    mid = [NaN NaN NaN];
    chanMidCache(key) = mid;
    return;
end

a = parts{1}; b = parts{2};

% Contacts already in canonical form like DRA1
if ~isKey(contactMap, a) || ~isKey(contactMap, b)
    mid = [NaN NaN NaN];
else
    %     xa = contactMap(a); xb = contactMap(b);
    %     mid = (xa + xb) / 2;
    xa = contactMap(a);   % already 1x3 numeric
    xb = contactMap(b);   % already 1x3 numeric
    mid = (xa + xb) / 2;


end

chanMidCache(key) = mid;
end

function pair_dist_mm = compute_pair_distance_mm(subj, pairNames_subject, coordsDir)

coordsPath = fullfile(coordsDir, subj + "_native_space_coords.xlsx");
if ~exist(coordsPath,"file")
    error("[DIST] Missing coords file for %s: %s", subj, coordsPath);
end

Tcoord = readtable(coordsPath);

% ---- Build monopolar contact key -> xyz ----
% IMPORTANT: THIS must match your FC contact naming (DLA1, DLH2, etc.)
% If your coord file stores contact as: ElectrodeName_labConvention_ + ContactNumber
% keep only rows that have coordinates
hasXYZ = ~isnan(Tcoord.x) & ~isnan(Tcoord.y) & ~isnan(Tcoord.z);
Tcoord = Tcoord(hasXYZ, :);

% ---- Build monopolar key -> xyz ----
base = upper(string(Tcoord.ElectrodeName_labConvention_));   % e.g. "dLA"
base = strrep(base,"_","");
base = strrep(base," ","");

% this makes "DLA" from "dLA" (since upper() already makes it "DLA")
% so we DON'T remove any letters; we just use the cleaned uppercase base.
contactKey = base + string(Tcoord.ContactNumber);            % e.g. "DLA1"

xyz = [Tcoord.x, Tcoord.y, Tcoord.z];

% build map
contactMap = containers.Map(cellstr(contactKey), num2cell(xyz,2));


% Cache midpoint per bipolar channel
chanMid = containers.Map('KeyType','char','ValueType','any');

pair_dist_mm = nan(numel(pairNames_subject),1);

for p = 1:numel(pairNames_subject)
    pid = char(pairNames_subject(p));  % "CHAN1--CHAN2"
    parts = strsplit(pid, "--");
    bip1 = parts{1}; bip2 = parts{2};

    m1 = get_midpoint_for_bipolar(bip1, contactMap, chanMid);
    m2 = get_midpoint_for_bipolar(bip2, contactMap, chanMid);

    if any(isnan(m1)) || any(isnan(m2))
        pair_dist_mm(p) = NaN;
        continue;
    end

    pair_dist_mm(p) = norm(m1 - m2);
end

% ============================
% NEW: Export bipolar node midpoints for Python plotting
% (uses chanMid cache you already built)
% ============================

% 1) collect unique bipolar node names from pairNames_subject ("BIP1--BIP2")
allBip = strings(0,1);
for p = 1:numel(pairNames_subject)
    parts = strsplit(char(pairNames_subject(p)), "--");
    if numel(parts) ~= 2, continue; end
    allBip(end+1,1) = string(parts{1});
    allBip(end+1,1) = string(parts{2});
end
allBip = unique(allBip);

% 2) compute/get midpoints (this will reuse cached values in chanMid)
midXYZ = nan(numel(allBip), 3);
for i = 1:numel(allBip)
    m = get_midpoint_for_bipolar(allBip(i), contactMap, chanMid);
    midXYZ(i,:) = m;
end

% 3) save to CSV
Tbip = table(allBip, midXYZ(:,1), midXYZ(:,2), midXYZ(:,3), ...
    'VariableNames', {'bipolar_name','x','y','z'});

outBipCsv = fullfile(coordsDir, subj + "_bipolar_nodes_midpoints.csv");
writetable(Tbip, outBipCsv);

fprintf("[BIPNODES] %s: saved %d bipolar midpoints to %s (NaN=%d)\n", ...
    subj, height(Tbip), outBipCsv, sum(any(isnan(midXYZ),2)));


fprintf("[DIST] %s: computed distances for %d edges (NaN=%d)\n", ...
    subj, numel(pair_dist_mm), sum(isnan(pair_dist_mm)));
end

