%% written by Tahereh Rashnavadi
% Jan 2026build_bipolar_signals_from_monopolar

function eegBinPath = find_eeg_bin_file(baseEEGDir, subjectID, runID)
%FIND_EEG_BIN_FILE Locate the IED-cleaned EEG .bin for a subject and run.
%
% baseEEGDir should be:
%   /Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis
%
% runID can be "Run2a" (fMRI segmentation) -> EEG typically stored under "Run2"
%
% Returns full path to:
%   <...>/<subj>/3_EEG/2_Cleaned/<subj>_RunX_Cleaned/IED_Cleaned/<subj>_RunX_IED_*Hz.bin

    runID = string(runID);
    baseRun = regexp(runID, 'Run\d+', 'match', 'once');  % "Run2" from "Run2a"
    if isempty(baseRun)
        error('find_eeg_bin_file:BadRunID', 'Cannot parse base run from runID = %s', runID);
    end
    baseRun = string(baseRun);

    eegDir = fullfile(baseEEGDir, subjectID, '3_EEG', '2_Cleaned', ...
        sprintf('%s_%s_Cleaned', subjectID, baseRun), 'IED_Cleaned');

    if ~isfolder(eegDir)
        error('find_eeg_bin_file:MissingDir', 'EEG directory not found: %s', eegDir);
    end

    % look for the IED cleaned bin
    patt = sprintf('%s_%s_IED_*Hz.bin', subjectID, baseRun);
    d = dir(fullfile(eegDir, patt));

    if isempty(d)
        error('find_eeg_bin_file:MissingBin', 'No EEG .bin found in %s matching %s', eegDir, patt);
    end

    % If multiple, pick the newest (or you can enforce 1)
    [~, idx] = max([d.datenum]);
    eegBinPath = fullfile(eegDir, d(idx).name);
end
