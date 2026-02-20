% Computes static EEG–fMRI similarity per band
% Finds the subject closest to the group-level profile
% Returns the statistically most representative subject

clear; clc;

% =========================
% PATHS
% =========================
fcRoot = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild";

bands = ["delta","theta","alpha","beta","gammalow","gammahigh"];
minPairs = 50;  % minimum number of edges required per band (adjust if needed)

% =========================
% FIND SUBJECTS
% =========================
D = dir(fcRoot);
isSubj = startsWith({D.name}, "ICE") & [D.isdir];
subjects = string({D(isSubj).name});
subjects = sort(subjects);

nSub  = numel(subjects);
nBand = numel(bands);

Rmat   = nan(nSub, nBand);   % Spearman similarity per subject x band
nEdge  = zeros(nSub, nBand); % number of edges used

fprintf("Found %d ICE folders in FC_rebuild.\n", nSub);

% =========================
% DEBUG SETTINGS
% =========================
DEBUG = true;
DEBUG_MAX_SUBJ = 10;     % full debug detail for first N subjects
NORM_NAMES = false;      % set true if name intersections are too small due to formatting
dbgShown = 0;

% IMPORTANT: Use EEG band field names that actually exist in your .mat files
bands = ["delta","theta","alpha","beta","gammalow","gammahigh"];
nBand = numel(bands);

% =========================
% MAIN LOOP
% =========================
for s = 1:nSub
    subj = subjects(s);

    eegFile  = fullfile(fcRoot, subj, "EEG_FC_Hilbert_" + subj + ".mat");
    fmriFile = fullfile(fcRoot, subj, "FMRI_FC_" + subj + ".mat");

    if ~isfile(eegFile) || ~isfile(fmriFile)
        if DEBUG
            fprintf("[SKIP] %s missing EEG or fMRI .mat\n", subj);
        end
        continue;
    end

    EEG  = load(eegFile,  "OUT");
    FMRI = load(fmriFile, "OUT");

    if ~isfield(EEG,"OUT") || ~isfield(FMRI,"OUT")
        if DEBUG
            fprintf("[SKIP] %s OUT missing in one file\n", subj);
        end
        continue;
    end

    % Pick a runKey that exists in BOTH EEG and fMRI
    runKey = pick_common_runkey(EEG.OUT, FMRI.OUT);
    if runKey == ""
        if DEBUG
            fprintf("[SKIP] %s no common runKey.\n", subj);
            fprintf("  EEG keys:  %s\n", strjoin(string(fieldnames(EEG.OUT))', ", "));
            fprintf("  fMRI keys: %s\n", strjoin(string(fieldnames(FMRI.OUT))', ", "));
        end
        continue;
    end

    % Show debug info for first few subjects only
    if DEBUG && dbgShown < DEBUG_MAX_SUBJ
        dbgShown = dbgShown + 1;
        fprintf("\n[DBG] %s using runKey=%s\n", subj, runKey);
        fprintf("  EEG.OUT.(runKey) fields: %s\n", strjoin(string(fieldnames(EEG.OUT.(runKey)))', ", "));
        fprintf("  fMRI.OUT.(runKey) fields: %s\n", strjoin(string(fieldnames(FMRI.OUT.(runKey)))', ", "));
    end

    % Get fMRI static z + names
    [Fz, fNames] = get_fmri_static_z(FMRI.OUT, runKey);
    if isempty(Fz) || isempty(fNames)
        if DEBUG
            fprintf("[FAIL] %s fMRI static matrix or names not found.\n", subj);
        end
        continue;
    end

    % Optional normalization if formatting differs
    if NORM_NAMES
        fNames = normalize_names(fNames);
    else
        fNames = string(fNames(:));
    end

    % Loop bands
    for b = 1:nBand
        band = bands(b);

        [Ez, eNames, why] = get_eeg_band_static_z_DEBUG(EEG.OUT, runKey, band);

        if isempty(Ez)
            if DEBUG && dbgShown <= DEBUG_MAX_SUBJ
                fprintf("  [MISS] %s band=%s EEG static missing: %s\n", subj, band, why);
            end
            continue;
        end

        if isempty(eNames)
            if DEBUG && dbgShown <= DEBUG_MAX_SUBJ
                fprintf("  [MISS] %s band=%s EEG names missing.\n", subj, band);
            end
            continue;
        end

        if NORM_NAMES
            eNames = normalize_names(eNames);
        else
            eNames = string(eNames(:));
        end

        % Align names
        [common, iF, iE] = intersect(fNames, eNames, 'stable');

        if numel(common) < 3
            if DEBUG && dbgShown <= DEBUG_MAX_SUBJ
                fprintf("  [MISS] %s band=%s name intersection too small: %d\n", subj, band, numel(common));
            end
            continue;
        end

        % Subset and vectorize
        Fzz = Fz(iF, iF);
        Ezz = Ez(iE, iE);

        fv = uppertri_vec(Fzz);
        ev = uppertri_vec(Ezz);

        good = ~isnan(fv) & ~isnan(ev) & isfinite(fv) & isfinite(ev);
        fv = fv(good);
        ev = ev(good);

        if numel(fv) < minPairs
            if DEBUG && dbgShown <= DEBUG_MAX_SUBJ
                fprintf("  [MISS] %s band=%s too few pairs after filtering: %d\n", subj, band, numel(fv));
            end
            continue;
        end

        Rmat(s,b)  = corr(fv, ev, 'Type','Spearman');
        nEdge(s,b) = numel(fv);
    end

    fprintf("[OK] %s (%s): bands computed = %d/%d\n", subj, runKey, sum(~isnan(Rmat(s,:))), nBand);
end

% =========================
% RANKING (ROBUST TO MISSING BANDS)
% =========================
minBandsRequired = 3;
valid = sum(~isnan(Rmat),2) >= minBandsRequired;

subjects_valid = subjects(valid);
R_valid = Rmat(valid,:);

fprintf("\nSubjects with >=%d bands available: %d\n", minBandsRequired, numel(subjects_valid));

if isempty(subjects_valid)
    warning("No subjects had enough bands computed. Likely name mismatch between EEG labels and fMRI seedNames. Try setting NORM_NAMES=true, or inspect one subject’s fNames/eNames.");
    return;
end

% Group profile (median per band across valid subjects)
groupProfile = nanmedian(R_valid, 1);

% Distance to group profile using available bands only
dist = nan(numel(subjects_valid),1);
for i = 1:numel(subjects_valid)
    mask = ~isnan(R_valid(i,:)) & ~isnan(groupProfile);
    if sum(mask) < minBandsRequired
        dist(i) = NaN;
    else
        dist(i) = sqrt(mean((R_valid(i,mask) - groupProfile(mask)).^2)); % normalized RMSE
    end
end

[dist_sorted, ord] = sort(dist, 'ascend', 'MissingPlacement','last');
rankedSubjects = subjects_valid(ord);

% --- ensure column vectors ---
rankedSubjects = rankedSubjects(:);
dist_sorted    = dist_sorted(:);

T = table(rankedSubjects, dist_sorted, ...
    'VariableNames', {'subject','dist_to_groupMedian'});

disp("=== Ranked subjects (most representative first) ===");
disp(T(1:min(20,height(T)),:));

fprintf("\nMost representative subject: %s\n", rankedSubjects(1));

% =========================
% -------- helpers --------
% =========================
function runKey = pick_common_runkey(OUTe, OUTf)
re = string(fieldnames(OUTe));
rf = string(fieldnames(OUTf));

re = re(startsWith(re,"Run"));
rf = rf(startsWith(rf,"Run"));

common = intersect(re, rf, 'stable');
if isempty(common)
    runKey = "";
else
    runKey = common(1); % choose the first common run (you can customize)
end
end

function [Fz, names] = get_fmri_static_z(OUT, runKey)
Fz = []; names = [];
if ~isfield(OUT, runKey); return; end
R = OUT.(runKey);

% names
if isfield(R,"seedNames")
    names = string(R.seedNames);
elseif isfield(R,"bipolarNames")
    names = string(R.bipolarNames);
else
    names = [];
end
names = names(:);

% matrix: prefer z if exists
if isfield(R,"FC_static_z")
    Fz = R.FC_static_z;
elseif isfield(R,"FC_static")
    % If only r exists, convert carefully (clip to avoid Inf)
    X = R.FC_static;
    X(1:size(X,1)+1:end) = NaN;        % ignore diag
    X = max(min(X, 0.999999), -0.999999);
    Fz = atanh(X);
    % put diag back to NaN
    Fz(1:size(Fz,1)+1:end) = NaN;
end
end

function [Ez, names] = get_eeg_band_static_z(OUT, runKey, band)
Ez = []; names = [];
if ~isfield(OUT, runKey); return; end
R = OUT.(runKey);

% names
if isfield(R,"seedNames")
    names = string(R.seedNames);
elseif isfield(R,"bipolarNames")
    names = string(R.bipolarNames);
else
    names = [];
end
names = names(:);

% band struct: OUT.Run1a.delta.FC_static_z style
bname = lower(string(band));
if isfield(R, bname)
    B = R.(bname);
    if isfield(B, "FC_static_z")
        Ez = B.FC_static_z;
        Ez(1:size(Ez,1)+1:end) = NaN; % ignore diag
        return;
    elseif isfield(B, "FC_static")
        X = B.FC_static;
        X(1:size(X,1)+1:end) = NaN;
        X = max(min(X, 0.999999), -0.999999);
        Ez = atanh(X);
        Ez(1:size(Ez,1)+1:end) = NaN;
        return;
    end
end
end

function v = uppertri_vec(M)
    n = size(M,1);
    mask = triu(true(n), 1);
    v = M(mask);
end

function [Ez, names, why] = get_eeg_band_static_z_DEBUG(OUT, runKey, band)
Ez = []; names = []; why = "";

if ~isfield(OUT, runKey)
    why = "OUT missing runKey";
    return;
end
R = OUT.(runKey);

% ---- names ----
if isfield(R,"seedNames")
    names = string(R.seedNames); names = names(:);
elseif isfield(R,"labels")
    % your EEG uses 'labels' (seen in debug output)
    names = string(R.labels); names = names(:);
elseif isfield(R,"bipolarNames")
    names = string(R.bipolarNames); names = names(:);
elseif isfield(R,"chanNames")
    names = string(R.chanNames); names = names(:);
else
    names = [];
end

bname = lower(string(band));

% ---- band struct: OUT.Run1a.delta, OUT.Run1a.gammalow, etc. ----
if ~isfield(R, bname)
    why = "Band field not found. (e.g., expected OUT.(runKey)."+bname+")";
    return;
end

B = R.(bname);

% ---- matrix fields inside band ----
if isfield(B, "FC_static_z")
    Ez = B.FC_static_z;
    Ez(1:size(Ez,1)+1:end) = NaN; % ignore diagonal
    return;
elseif isfield(B, "FC_static")
    X = B.FC_static;
    X(1:size(X,1)+1:end) = NaN;
    X = max(min(X, 0.999999), -0.999999);
    Ez = atanh(X);
    Ez(1:size(Ez,1)+1:end) = NaN;
    return;
else
    why = "Band struct exists but missing FC_static_z/FC_static fields.";
    return;
end
end

