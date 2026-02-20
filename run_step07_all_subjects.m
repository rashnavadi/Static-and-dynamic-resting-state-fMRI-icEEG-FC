% run_step07_all_subjects.m

clear; clc;

% ---- Paths ----
projectRoot = '/Volumes/MIND/ICE/Tara/Kristina_paper';
fcRoot      = fullfile(projectRoot, 'FC_rebuild');

masterOutDir = fullfile(projectRoot, 'master_tables');
demoCsvPath  = fullfile(projectRoot, 'demographics.csv');

makeWindowTable = true;   % or false if you only want pair-level


% ---- Collect all config files ----
cfgFiles = dir(fullfile(fcRoot, 'ICE*', 'config_ICE*.mat'));

% ---- Keep only a subject range (e.g., ICE050 to ICE070) ----
minID = 1;
maxID = 71;

keep = false(size(cfgFiles));

for i = 1:numel(cfgFiles)
    fname = cfgFiles(i).name;  % e.g., "config_ICE062.mat"
    tok = regexp(fname, 'ICE(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        idNum = str2double(tok{1});
        keep(i) = (idNum >= minID) && (idNum <= maxID);
    end
end

cfgFiles = cfgFiles(keep);

fprintf("Keeping %d config files (ICE%03d–ICE%03d)\n", numel(cfgFiles), minID, maxID);

cfgPaths = cell(numel(cfgFiles), 1);
for i = 1:numel(cfgFiles)
    cfgPaths{i} = fullfile(cfgFiles(i).folder, cfgFiles(i).name);
end

fprintf('Found %d subjects\n', numel(cfgPaths));

% ---- Run Step 07 ----
Step07_build_master_tables(cfgPaths, masterOutDir, demoCsvPath, makeWindowTable);
