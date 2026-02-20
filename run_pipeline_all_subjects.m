%% Run ICE pipeline for all subjects
% Written by Tahereh Rashnavadi
% Jan 2026

clear; clc;

% ---- Paths ----
projectRoot = "/Volumes/MIND/ICE/Tara/Kristina_paper";
fmriTSRoot  = fullfile(projectRoot, "fMRI_timeseries");

pipelineDir = "/Volumes/MIND/ICE/scripts/Kristina_paper";
addpath(pipelineDir);
addpath(fullfile(pipelineDir, "helpers"));

% ---- Get subject list automatically ----
d = dir(fmriTSRoot);
d = d([d.isdir]);
subjects = string({d.name});
subjects = subjects(startsWith(subjects, "ICE"));

fprintf("Found %d subjects\n", numel(subjects));
disp(subjects');

% ---- Loop over subjects ----
for s = 1:numel(subjects)

    subj = subjects(s);
    subj_safe = subj;   % <-- PROTECT FROM clear

    fprintf("\n==============================\n");
    fprintf("Running pipeline for %s\n", subj_safe);
    fprintf("==============================\n");

    try
        % ---------- STEP 00 ----------
        C = struct();
        C.subjectID = subj_safe;

        run("step00_config.m");

%         % ---------- STEP 01–05 ----------
%         run("step01_make_segment_table.m");
%         run("step02_load_fmri_seedTS.m");
%         run("step03_load_eeg_and_segment.m");
%         run("step04_compute_fc_fmri.m");
        run("step05_compute_fc_eeg_hilbert_use_fmri_final_seeds.m");

        fprintf("[DONE] %s completed successfully\n", subj_safe);

    catch ME
        fprintf(2, "[ERROR] %s failed\n", subj_safe);
        fprintf(2, "%s\n", ME.message);
        continue;
    end
end

fprintf("\nALL SUBJECTS FINISHED\n");


