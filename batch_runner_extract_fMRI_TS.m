%% run_extract_bipolar_fMRI_seedTS_excludeIED_ALL.m
% Runs extract_bipolar_fMRI_seedTS_excludeIED for all subjects and all runs
% found under: /Volumes/MIND/ICE/Tara/fMRI_timeseries/<SUBJECT>/
%
% Tahereh Rashnavadi - Jan 2026

clear; clc;

% ---- Fixed paths (same as your single-subject call) ----
baseFMriDir   = "/Volumes/MIND/ICE/ICE_denoised_filtered_funcs";
baseElecDir   = "/Volumes/MIND/ICE/Tara/native_space_electrodes_coords";
mainMapMat    = "/Volumes/MIND/ICE/scripts/Kristina_paper/helpers/main_channels_map_bipolar.mat";
projectRoot   = "/Volumes/MIND/ICE/Tara/Kristina_paper";

% Where runs are discovered from (folders like Run1a, Run2a, Run3, ...)
baseRunListDir = "/Volumes/MIND/ICE/Tara/fMRI_timeseries";

% ---- Subject list ----
% subjects = [ ...
% "ICE013","ICE017","ICE023","ICE028","ICE031","ICE035","ICE038","ICE041","ICE044","ICE047","ICE050","ICE053","ICE056","ICE059","ICE063","ICE066", ...
% "ICE014","ICE018","ICE024","ICE029","ICE033","ICE036","ICE039","ICE042","ICE045","ICE048","ICE051","ICE054","ICE057","ICE060","ICE064","ICE069", ...
% "ICE016","ICE022","ICE027","ICE030","ICE034","ICE037","ICE040","ICE043","ICE046","ICE049","ICE052","ICE055","ICE058","ICE062","ICE065","ICE070" ...
% ];
subjects = [ ...
"ICE023","ICE028","ICE031","ICE035","ICE038","ICE041","ICE044","ICE047","ICE050","ICE053","ICE056","ICE059","ICE063","ICE066", ...
"ICE024","ICE029","ICE033","ICE036","ICE039","ICE042","ICE045","ICE048","ICE051","ICE054","ICE057","ICE060","ICE064","ICE069", ...
"ICE016","ICE022","ICE027","ICE030","ICE034","ICE037","ICE040","ICE043","ICE046","ICE049","ICE052","ICE055","ICE058","ICE062","ICE065","ICE070" ...
];

% ---- Optional: add your code/helpers to path if needed ----
% addpath("/Volumes/MIND/ICE/scripts/Kristina_paper");
% addpath(genpath("/Volumes/MIND/ICE/scripts/Kristina_paper/helpers"));

% ---- Logging ----
logDir = fullfile(projectRoot, "batch_logs");
if ~exist(logDir, "dir"), mkdir(logDir); end
logFile = fullfile(logDir, sprintf("batch_extract_seedTS_%s.txt", datestr(now,"yyyymmdd_HHMMSS")));
fid = fopen(logFile, "w");
fprintf(fid, "Batch run started: %s\n", datestr(now));
fprintf(fid, "baseRunListDir: %s\n\n", baseRunListDir);

% Keep a structured summary in memory too
summary = table('Size',[0 5], ...
    'VariableTypes',{'string','string','string','string','datetime'}, ...
    'VariableNames',{'Subject','Run','Status','Message','Timestamp'});

% ---- Main loop ----
for s = 1:numel(subjects)
    subjectID = subjects(s);
    subjRunDir = fullfile(baseRunListDir, subjectID);

    if ~isfolder(subjRunDir)
        msg = "Run directory not found (skipping subject).";
        fprintf("❌ %s: %s\n", subjectID, msg);
        fprintf(fid, "❌ %s: %s | %s\n", subjectID, msg, subjRunDir);
        summary = [summary; {subjectID,"", "SKIP_SUBJECT", msg, datetime('now')}]; %#ok<AGROW>
        continue;
    end

    d = dir(subjRunDir);
    runNames = string({d([d.isdir]).name});
    runNames = runNames(~ismember(runNames, [".",".."]));

    % Only keep folders that look like Run*
    isRun = startsWith(runNames, "Run");
    runNames = runNames(isRun);

    % Natural-ish sorting: Run1, Run1a, Run2, Run10, ...
    runNames = sort_run_ids(runNames);

    if isempty(runNames)
        msg = "No Run* folders found (skipping subject).";
        fprintf("⚠️  %s: %s\n", subjectID, msg);
        fprintf(fid, "⚠️  %s: %s | %s\n", subjectID, msg, subjRunDir);
        summary = [summary; {subjectID,"", "SKIP_SUBJECT", msg, datetime('now')}]; %#ok<AGROW>
        continue;
    end

    fprintf("\n===== %s (%d runs) =====\n", subjectID, numel(runNames));
    fprintf(fid, "\n===== %s (%d runs) =====\n", subjectID, numel(runNames));

    for r = 1:numel(runNames)
        runID = runNames(r);

        try
            fprintf("➡️  Running %s %s ...\n", subjectID, runID);
            fprintf(fid, "➡️  Running %s %s ...\n", subjectID, runID);

            extract_bipolar_fMRI_seedTS_excludeIED( ...
                char(subjectID), char(runID), ...
                char(baseFMriDir), ...
                char(baseElecDir), ...
                char(mainMapMat), ...
                char(projectRoot) );

            fprintf("✅ Done %s %s\n", subjectID, runID);
            fprintf(fid, "✅ Done %s %s\n", subjectID, runID);

            summary = [summary; {subjectID, runID, "OK", "Completed", datetime('now')}]; %#ok<AGROW>

        catch ME
            errMsg = string(ME.message);
            fprintf("❌ FAILED %s %s | %s\n", subjectID, runID, errMsg);
            fprintf(fid, "❌ FAILED %s %s | %s\n", subjectID, runID, errMsg);

            % also log stack (first few lines) for debugging
            st = ME.stack;
            for k = 1:min(5, numel(st))
                fprintf(fid, "    at %s (line %d)\n", st(k).name, st(k).line);
            end

            summary = [summary; {subjectID, runID, "FAIL", errMsg, datetime('now')}]; %#ok<AGROW>

            % Keep going to next run/subject
        end
    end
end

fprintf(fid, "\nBatch run finished: %s\n", datestr(now));
fclose(fid);

% Save summary table
summaryFile = fullfile(logDir, sprintf("batch_extract_seedTS_summary_%s.csv", datestr(now,"yyyymmdd_HHMMSS")));
writetable(summary, summaryFile);

fprintf("\n✅ Batch finished.\nLog: %s\nSummary: %s\n", logFile, summaryFile);

%% -------- Local helper: sort Run IDs like Run1, Run1a, Run2, Run10, Run10b ----------
function runNamesSorted = sort_run_ids(runNames)
% runNames: string array like ["Run1a","Run10","Run2","Run3","Run1"]
% Sorts by numeric part then by suffix ("" < "a" < "b" ...)

n = numel(runNames);
runNum = nan(n,1);
runSuf = strings(n,1);

for i = 1:n
    tok = regexp(runNames(i), '^Run(\d+)([A-Za-z]*)$', 'tokens', 'once');
    if isempty(tok)
        runNum(i) = inf;
        runSuf(i) = runNames(i);
    else
        runNum(i) = str2double(tok{1});
        runSuf(i) = lower(string(tok{2}));
    end
end

T = table(runNames(:), runNum, runSuf, 'VariableNames', {'name','num','suf'});
T = sortrows(T, {'num','suf','name'});
runNamesSorted = T.name;
end
