% Run it (example):
% export_FMRI_run_to_csv("ICE070","Run1");

function export_FMRI_run_to_csv(subjID, runKey)
% Exports fMRI seed names + FC_static (prefers z if available) as labeled CSV.

rootDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild";
subjDir = fullfile(rootDir, subjID);
fmriFile = fullfile(subjDir, "FMRI_FC_" + subjID + ".mat");
load(fmriFile, "OUT");

R = OUT.(runKey);

% ---- names ----
if isfield(R, "seedNames")
    names = string(R.seedNames);
elseif isfield(R, "labels")
    names = string(R.labels);
else
    error("Cannot find seed/channel names in OUT.%s (expected seedNames or labels).", runKey);
end

% ---- matrix (prefer Fisher z if you have it; otherwise raw r) ----
FC = pickField(R, ["FC_static_z","FC_static","FC_static_r","R_static_z","R_static","r_static"]);

% ---- output folder ----
outDir = fullfile(rootDir, "RESULTS_REPRESENTATIVE", "python_inputs", subjID);
if ~exist(outDir,"dir"), mkdir(outDir); end

% 1) seed names
outNames = fullfile(outDir, subjID + "_" + runKey + "_seedNames.csv");
writetable(table(names, 'VariableNames', {'seedNames'}), outNames);

% 2) labeled FC
T = array2table(FC, 'RowNames', cellstr(names), 'VariableNames', cellstr(names));
outFC = fullfile(outDir, subjID + "_" + runKey + "_FMRI_FC_static_labeled.csv");
writetable(T, outFC, 'WriteRowNames', true);

fprintf("[OK] fMRI names: %s\n", outNames);
fprintf("[OK] fMRI FC:    %s\n", outFC);
fprintf("[INFO] N=%d\n", numel(names));
end

function A = pickField(S, candidates)
for k = 1:numel(candidates)
    f = candidates(k);
    if isfield(S, f)
        A = S.(f);
        return;
    end
end
error("None of the FC fields found. Tried: %s", strjoin(cellstr(candidates), ", "));
end
