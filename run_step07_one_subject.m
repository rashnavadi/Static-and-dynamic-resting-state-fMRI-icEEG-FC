% run_step07_all_subjects.m
clear; clc;

% ---- Paths ----
projectRoot = '/Volumes/MIND/ICE/Tara/Kristina_paper';
fcRoot      = fullfile(projectRoot, 'FC_rebuild');

masterOutDir = fullfile(projectRoot, 'master_tables');
demoCsvPath  = fullfile(projectRoot, 'demographics.csv');

makeWindowTable = true;

% ---- Find ALL config files ----
cfgFiles = dir(fullfile(fcRoot, 'ICE*', 'config_ICE*.mat'));

if isempty(cfgFiles)
    error("No config_ICE*.mat files found under: %s", fcRoot);
end

% ---- Sort by ICE number (so it runs ICE001, ICE002, ...) ----
ids = nan(numel(cfgFiles), 1);
for i = 1:numel(cfgFiles)
    tok = regexp(cfgFiles(i).name, 'ICE(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        ids(i) = str2double(tok{1});
    end
end
[~, idx] = sort(ids);
cfgFiles = cfgFiles(idx);

% ---- Build full paths ----
cfgPaths = cell(numel(cfgFiles), 1);
for i = 1:numel(cfgFiles)
    cfgPaths{i} = fullfile(cfgFiles(i).folder, cfgFiles(i).name);
end

fprintf('Running Step07 for %d subjects\n', numel(cfgPaths));

% ---- Run Step07 ----
Step07_build_master_tables(cfgPaths, masterOutDir, demoCsvPath, makeWindowTable);
