% pipeline_ICE033/
%   00_config_ICE033.m
%   01_make_segment_table_ICE033.m
%   02_load_fmri_seedTS_ICE033.m
%   03_load_eeg_and_segment_ICE033.m
%   04_compute_fc_fmri_ICE033.m
%   05_compute_fc_eeg_hilbert_ICE033.m
%   06_export_for_R_ICE033.m
% 
%   helpers/
%     parse_segments_row.m
%     find_eeg_bin_file.m
%     load_eeg_bin_with_labels.m       
%     canonical_bipolar_name.m
%     match_bipolar_labels.m
%     bandpass_hilbert_envelope.m
%     uppertri_vec.m
%     sliding_window_indices.m
%% written by Tahereh Rashnavadi
% Jan 2026

% clear; clc;

% C = struct();

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);


% ---- code dir ----
C.codeDir = "/Volumes/MIND/ICE/scripts/Kristina_paper";  % where pipeline scripts + helpers live

addpath(fullfile(C.codeDir, "helpers"));

%% step00_config_ICE033.m

% ---- Subject / TR ----
if ~isfield(C, "subjectID") || isempty(C.subjectID)
    error("C.subjectID must be set by the calling script");
end
C.TR = 1.5;                 % verified constant
C.fmriLagSec = 6.0;         % Pierre/old code: use 6s so it aligns to TR
C.winSec = 120;
C.overlapSec = 60;

% ---- Base paths ----
C.baseFMriDir = "/Volumes/MIND/ICE/ICE_denoised_filtered_funcs";  % contains ICE033/Run2a/...
C.baseEEGDir  = "/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis"; % contains ICE033/3_EEG/...

% ---- Project root (where your extracted TS lives) ----
C.projectRoot = "/Volumes/MIND/ICE/Tara/Kristina_paper";

% ---- segments ----
C.segmentsMat = fullfile(C.codeDir, "helpers", "segments.mat");

% ---- fMRI time series (pre-extracted) ----
% Update these 2 patterns to your real storage location:
C.fmriTSDir = fullfile(C.projectRoot, "fMRI_timeseries", C.subjectID);
C.bipTS_pattern  = "%s_%s_seedTS_bipolar_exclIED.mat";

% ---- fMRI Runs (AUTO from folders inside fmriTSDir) ----
if ~exist(C.fmriTSDir, 'dir')
    error("Subject folder not found: %s", C.fmriTSDir);
end

d = dir(fullfile(C.fmriTSDir, "Run*"));
d = d([d.isdir]);                 % keep only directories
names = string({d.name})';

% remove '.' and '..' just in case, and keep only things that start with "Run"
names = names(~ismember(names, [".",".."]));
names = names(startsWith(names, "Run"));

if isempty(names)
    error("No run folders found in %s (expected Run*)", C.fmriTSDir);
end

% Natural-ish sort: Run1, Run1a, Run1b, Run2, ...
[~, idx] = sort(cellfun(@run_sort_key, cellstr(names)));
C.runs = names(idx);

fprintf("[OK] Auto-detected runs for %s: %s\n", C.subjectID, strjoin(C.runs, ", "));

% ---- Output ----
C.outDir = fullfile(C.projectRoot, "FC_rebuild", C.subjectID);
if ~exist(C.outDir, 'dir'), mkdir(C.outDir); end

% ---- LMM output (statistics) ----
C.lmmOutDir = fullfile(C.projectRoot, "lmm_output", C.subjectID);
if ~exist(C.lmmOutDir, 'dir'), mkdir(C.lmmOutDir); end


save(fullfile(C.outDir, sprintf("config_%s.mat", C.subjectID)), "C");
fprintf("[OK] Config saved: %s\n", fullfile(C.outDir, sprintf("config_%s.mat", C.subjectID)));



% ---------- helper for sorting ----------
function k = run_sort_key(r)
% r like 'Run1', 'Run5a', 'Run10b'
    m = regexp(r, '^Run(\d+)([a-z]?)$', 'tokens', 'once');
    if isempty(m)
        k = 1e9;  % push unknown names to end
        return;
    end
    num = str2double(m{1});
    suf = m{2};
    if isempty(suf)
        sufVal = 0;
    else
        sufVal = double(suf) - double('a') + 1; % a=1, b=2, ...
    end
    k = num*10 + sufVal;
end
