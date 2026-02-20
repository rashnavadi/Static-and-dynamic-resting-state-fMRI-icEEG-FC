
% HOW TO RUN IT:
% partsDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/parts_window";
% outDir   = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial";
% step09_spatialStratified_windowLMM(partsDir, outDir);

function step09_spatialStratified_windowLMM(partsDir, outDir)


if nargin < 1 || isempty(partsDir)
    partsDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/parts_window";
end
if nargin < 2 || isempty(outDir)
    outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial";
end
if ~exist(outDir,"dir"); mkdir(outDir); end

bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

% ---- choose bins (8 bins) ----
regions = ["LF","RF","LMT","RMT","LLT","RLT","LI","RI"];
nR = numel(regions);

% ---------- Load + concat ----------
files = dir(fullfile(partsDir, "ICE*_forLMM_window.mat"));
assert(~isempty(files), "No *_forLMM_window.mat found in %s", partsDir);

Tall = table();
for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name), "Tsub");
    Tall = [Tall; S.Tsub]; %#ok<AGROW>
end

% Put this EXACTLY here if you want it:
disp(Tall.Properties.VariableNames);

% ---------- Basic type fixes ----------
Tall.subject_id = categorical(string(Tall.subject_id));
Tall.pair_id    = string(Tall.pair_id);
Tall.band       = categorical(string(Tall.band), cellstr(bandOrder), "Ordinal", true);
Tall.run_id     = categorical(string(Tall.run_id));
Tall.sex        = categorical(string(Tall.sex));

% ---------- Parse pair_id => nodeA/nodeB ----------
[nodeA, nodeB] = split_pair_id_into_nodes(Tall.pair_id);
Tall.nodeA = nodeA;
Tall.nodeB = nodeB;
fprintf("\n[DEBUG] Example parsed pair_id → nodeA / nodeB:\n");
disp(Tall(1:5, {'pair_id','nodeA','nodeB'}));

% ---------- Map each node => region bin ----------
Tall.regA = map_node_to_region(Tall.nodeA);
Tall.regB = map_node_to_region(Tall.nodeB);

fprintf("\n[DEBUG] Unique region labels found:\n");
disp(categories(categorical([Tall.regA; Tall.regB]))');

fprintf("[DEBUG] Example node → region mapping:\n");
disp(Tall(1:8, {'nodeA','regA','nodeB','regB'}));

% Identify rows we can’t classify
bad = ismissing(Tall.regA) | ismissing(Tall.regB) | Tall.regA=="UNK" | Tall.regB=="UNK";

% DEBUG UNK nodes BEFORE dropping them
if any(bad)
    fprintf("[INFO] Found %d rows with UNK region labels (will drop)\n", sum(bad));
    unkNodes = unique([Tall.nodeA(bad); Tall.nodeB(bad)]);
    fprintf("[DEBUG] Example UNK nodes (first 30 unique):\n");
    disp(unkNodes(1:min(30,numel(unkNodes))));
end

% Now drop them
if any(bad)
    Tall(bad,:) = [];
end

% Show which region labels remain
fprintf("[DEBUG] Remaining region labels after dropping UNK:\n");
disp(categories(categorical([Tall.regA; Tall.regB]))');









% ---------- Ensure distance covariate exists ----------
needDist = ismember("c_log_pair_dist", Tall.Properties.VariableNames);
if ~needDist
    % fallback: if only pair_dist_mm exists
    assert(ismember("pair_dist_mm", Tall.Properties.VariableNames), "No distance column found in Tall.");
    Tall.pair_dist_mm = double(Tall.pair_dist_mm);
    Tall.log_pair_dist = log(Tall.pair_dist_mm);
    Tall.c_log_pair_dist = Tall.log_pair_dist - mean(Tall.log_pair_dist, 'omitnan');
end

% ---------- Model (same as your main), just add distance covariate ----------
formula = ...
 "fmri_fc_windowed ~ eeg_fc_trueStatic + delta_eeg_fc + eeg_fc_windowed_sd + c_log_pair_dist + age + sex" + ...
 " + (1|subject_id) + (1|subject_id:run_id) + (1|pair_id)";

% Thresholds so you don’t fit nonsense subsets
minObs   = 2000;
minSubj  = 5;
minPairs = 50;

for b = bandOrder
    fprintf("\n==============================\n");
    fprintf("[DEBUG] Processing band: %s\n", string(b));

    Tb0 = Tall(Tall.band == b, :);
    Tb0.subject_id = removecats(Tb0.subject_id);
    Tb0.run_id     = removecats(Tb0.run_id);
    Tb0.sex        = removecats(Tb0.sex);


    % output matrices for this band
    BETA_dev = nan(nR,nR);
    P_dev    = nan(nR,nR);
    BETA_static = nan(nR,nR);
    P_static    = nan(nR,nR);
    Nobs     = zeros(nR,nR);
    Nsubj    = zeros(nR,nR);
    Npairs   = zeros(nR,nR);

    for i = 1:nR
        for j = 1:nR
            A = regions(i); B = regions(j);

            % symmetric selection (LF–RF == RF–LF)
            sel = (Tb0.regA==A & Tb0.regB==B) | (Tb0.regA==B & Tb0.regB==A);
            Tb = Tb0(sel,:);

            % required columns
            useCols = ["fmri_fc_windowed","eeg_fc_trueStatic","delta_eeg_fc","eeg_fc_windowed_sd", ...
                       "c_log_pair_dist","age","sex","subject_id","run_id","pair_id"];
            Tb = rmmissing(Tb, "DataVariables", useCols);

            nObs = height(Tb);
            nS = numel(unique(string(Tb.subject_id)));
            nP = numel(unique(string(Tb.pair_id))); % unique pairs in this subset

            fprintf("\n[DEBUG] Band %s | %s–%s\n", string(b), A, B);
            fprintf("        nObs=%d | nSubj=%d | nPairs=%d\n", ...
                height(Tb), ...
                numel(unique(string(Tb.subject_id))), ...
                numel(unique(string(Tb.pair_id))));

            Nobs(i,j)=nObs; Nsubj(i,j)=nS; Npairs(i,j)=nP;

            if nObs < minObs || nS < minSubj || nP < minPairs
                continue;
            end

            try
                lme = fitlme(Tb, formula, "FitMethod","REML");
                r = lmm_extract(lme, "delta_eeg_fc");   % dynamic coupling term
                BETA_dev(i,j) = r.beta;
                P_dev(i,j)    = r.p;

                % static fc
                rS = lmm_extract(lme, "eeg_fc_trueStatic");
                BETA_static(i,j) = rS.beta;
                P_static(i,j)    = rS.p;

            catch ME
                fprintf("[WARN] Band %s %s-%s failed: %s\n", string(b), A, B, ME.message);
            end
        end
    end

    % Save per-band tables (matrix as table with row/col labels)
    Tbeta = array2table(BETA_dev, "VariableNames", cellstr(regions), "RowNames", cellstr(regions));
    Tp    = array2table(P_dev,    "VariableNames", cellstr(regions), "RowNames", cellstr(regions));
    Tnobs = array2table(Nobs,     "VariableNames", cellstr(regions), "RowNames", cellstr(regions));
    Tnsub = array2table(Nsubj,    "VariableNames", cellstr(regions), "RowNames", cellstr(regions));
    Tnp   = array2table(Npairs,   "VariableNames", cellstr(regions), "RowNames", cellstr(regions));
    TbetaS = array2table(BETA_static, "VariableNames", cellstr(regions), "RowNames", cellstr(regions));
    TpS    = array2table(P_static,    "VariableNames", cellstr(regions), "RowNames", cellstr(regions));

    writetable(Tbeta, fullfile(outDir, "BETA_dev_" + string(b) + ".csv"), "WriteRowNames", true);
    writetable(Tp,    fullfile(outDir, "P_dev_"    + string(b) + ".csv"), "WriteRowNames", true);
    writetable(Tnobs, fullfile(outDir, "Nobs_"     + string(b) + ".csv"), "WriteRowNames", true);
    writetable(Tnsub, fullfile(outDir, "Nsubj_"    + string(b) + ".csv"), "WriteRowNames", true);
    writetable(Tnp,   fullfile(outDir, "Npairs_"   + string(b) + ".csv"), "WriteRowNames", true);
    writetable(TbetaS, fullfile(outDir, "BETA_static_" + string(b) + ".csv"), "WriteRowNames", true);
    writetable(TpS,    fullfile(outDir, "P_static_"    + string(b) + ".csv"), "WriteRowNames", true);

    fprintf("[OK] Saved spatial matrices for band %s\n", string(b));
end

fprintf("[DONE] Spatially stratified window LMM outputs in:\n%s\n", outDir);
end

function [nodeA, nodeB] = split_pair_id_into_nodes(pair_id)
pair_id = string(pair_id);
parts = split(pair_id, "--");
nodeA = parts(:,1);
nodeB = parts(:,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function reg = map_node_to_region(nodeName)

nodeName = string(nodeName);

% nodeName is bipolar like "DRAIN1-DRAIN2"
% we convert to a single "base label" like "DRAIN" using the FIRST contact
base = extractBefore(nodeName, "-");   % "DRAIN1"
base = regexprep(base, "\d+$", "");    % remove trailing digits -> "DRAIN"
base = upper(base);

reg = strings(size(base));

for k = 1:numel(base)
    b = base(k);

    % laterality
    if contains(b, "L")
        lat = "L";
    elseif contains(b, "R")
        lat = "R";
    else
        reg(k) = "UNK";
        continue;
    end

    % Insula: IN anywhere
    if contains(b, "IN")
        reg(k) = lat + "I";    % LI / RI
        continue;
    end

    % Mesial temporal: hippocampus/amygdala/uncus
    % (adjust these if your codes differ)
    if contains(b, "H") || contains(b, "A") || contains(b, "U")
        reg(k) = lat + "MT";   % LMT / RMT
        continue;
    end

    % Temporal (but not mesial)
    if contains(b, "T") || contains(b, "HG")
        reg(k) = lat + "LT";   % LLT / RLT
        continue;
    end

    % Frontal: F/PF/OF/OP/SMA/M1/etc.
    if contains(b, "F") || contains(b, "PF") || contains(b, "OF") || contains(b, "OP") || contains(b, "SMA") || contains(b, "M1")
        reg(k) = lat + "F";    % LF / RF
        continue;
    end

    % If you want to add Parietal/Occipital later:
    % if contains(b,"P") => LP/RP
    % if contains(b,"O") => LO/RO

    reg(k) = "UNK";
end

reg = categorical(reg);
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

