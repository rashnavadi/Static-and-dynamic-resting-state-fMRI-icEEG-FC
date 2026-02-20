function main_channels_map_bipolar = main_channels_map(baseDataDir)
% main_bipolar_channels_map
%
% Builds a subject/IED-specific map of *bipolar* main channels based on:
%   1) A hard-coded list of MONOPOLAR main contacts per subject/IED
%   2) The actual existing monopolar channels from each subject's
%      ProcessingLog_&_ChannelLabels.txt (Run1)
%   3) Generating only bipolar neighbors that truly exist
%
% Keys look like: 'ICE017_IED1'
% Values are cell arrays of bipolar labels: {'DRAIN2-DRAIN3', 'DRAIN3-DRAIN4', ...}
%
% Usage:
%   main_map = main_bipolar_channels_map();
%   main_map = main_bipolar_channels_map('/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis');
%
% The function also saves:
%   main_channels_map_bipolar.mat
%
% Tara / Pierre ROC/FC project

    if nargin < 1
        baseDataDir = '/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis';
    end

    % ---------- 1. Define MONOPOLAR main channels (per subject/IED) ----------
    mono_main = containers.Map();

    mono_main('ICE013_IED1') = {'dLH4', 'dLH5'};
    % mono_main('ICE013_IED2') = {'gLiT5'};
    % mono_main('ICE013_IED3') = {'sLpT1', 'sLpT2', 'sLpT3', 'sLpT4'};

    mono_main('ICE014_IED1') = {'dRaIN3', 'dRaIN5'}; % negative
    mono_main('ICE014_IED2') = {'dRmIN3', 'dRmIN4'}; % positive
    mono_main('ICE014_IED3') = {'dLpIN1', 'dLpIN2'}; % positive
    mono_main('ICE014_IED4') = {'dLmT2', 'dLmT3'};   % positive
    mono_main('ICE014_IED5') = {'dRmT5', 'dRmT6'};   % positive

    % ICE015 excluded

    mono_main('ICE016_IED1') = {'dLmH1'}; % dLaH1 removed due to noise

    mono_main('ICE017_IED1') = {'dRpT1', 'dRpT2'};   % negative
    mono_main('ICE017_IED2') = {'dRH3', 'dRH4'};     % negative

    mono_main('ICE018_IED1') = {'dRH1', 'dRH2'};
    mono_main('ICE018_IED2') = {'dLH1', 'dLH2'};

    % ICE019, ICE020, ICE021 custom/excluded in your notes – omitted here

    mono_main('ICE022_IED1') = {'dLA3', 'dLA4'};    % variable polarity
    mono_main('ICE022_IED2') = {'dRH1', 'dRH2'};    % variable polarity

    mono_main('ICE023_IED1') = {'dRmsT5', 'dRmsT6', 'dRpsT5', 'dRpsT6'};

    mono_main('ICE024_IED1') = {'dLpIN1', 'dLpIN2', 'dLpIN3'};

    mono_main('ICE027_IED1') = {'dLasT4', 'dLmsT2', 'dLmsT3'};
    mono_main('ICE027_IED2') = {'dLA7', 'dLA8'};

    mono_main('ICE028_IED1') = {'dLaR1', 'dLpR4', 'dLpR5', 'dLpR6'};
    mono_main('ICE028_IED2') = {'dLaR1', 'dLpR4', 'dLpR5', 'dLpR6'};

    mono_main('ICE029_IED1') = {'dRaH1', 'dRaH2', 'dRaH3', 'dRaH4', 'dRpH1', 'dRpH2', 'dRpH3'};
    mono_main('ICE029_IED2') = {'dLaH1', 'dLaH2', 'dLaH3', 'dLA2', 'dLA3', 'dLA4'};

    mono_main('ICE030_IED1') = {'dLaH1', 'dLaH7', 'dLaH8', 'dLpH1', 'dLpH7', 'dLpH8'};

    mono_main('ICE031_IED1') = {'dLSMA6', 'dLSMA7', 'dLSMA8'};
    mono_main('ICE031_IED2') = {'dLPM5', 'dLPM6', 'dLPM7'};

    mono_main('ICE033_IED1') = {'dRasTg1', 'dRasTg2', 'dRasTg3', 'dRasTg4', 'dRaIN3', 'dRaIN4'};
    mono_main('ICE033_IED2') = {'dRlOF2', 'dRlOF3', 'dRlOF4', 'dRlOF5', 'dRlOF7', 'dRaIN3', 'dRaIN4'};

    mono_main('ICE034_IED1') = {'dLA1', 'dLA2', 'dLaH1', 'dLaH2', 'dLaH3'};
    mono_main('ICE034_IED2') = {'dRaH1', 'dRaH2', 'dRA2', 'dRA3'};

    mono_main('ICE035_IED1') = {'dRmF3', 'dRmF4', 'dRH5', 'dRH6', 'dRaIN1', 'dRaIN2', 'dRaIN3', 'dRaIN4'};

    mono_main('ICE036_IED1') = {'dRpIN1', 'dRpIN2', 'dRpIN3', 'dRpIN4', 'dRpIN5', ...
                                 'dRaIN1', 'dRaIN2', 'dRaIN3', 'dRaIN4', 'dRaIN5', ...
                                 'dRH1', 'dRH2', 'dRH3'};
    mono_main('ICE036_IED2') = {'dLH1', 'dLH2', 'dLH3'};

    mono_main('ICE037_IED1') = {'dRA2', 'dRA3', 'dRA4', 'dRaH1', 'dRaH2', 'dRaH3'};

    mono_main('ICE038_IED1') = {'dLaH3', 'dLaH4', 'dLaH5', 'dLaH6', 'dLpH1', 'dLpH2', 'dLpH3', 'dLpH4'};
    mono_main('ICE038_IED2') = {'dLmO7', 'dLmO8', 'dLmO9', 'dLmO10'};

    mono_main('ICE039_IED1') = {'dLaxH3', 'dLaxH4', 'dLaxH5', 'dLaxH6'};
    mono_main('ICE039_IED2') = {'dLamTg7', 'dLamTg8', 'dLamTg9', ...
                                 'dLmmTg7', 'dLmmTg8', 'dLmmTg9', 'dLmmTg10', ...
                                 'dLpmTg5', 'dLpmTg6', 'dLpmTg7', 'dLpmTg8'};

    mono_main('ICE040_IED1') = {'dRaH3', 'dRaH4', 'dRaH5', 'dRaH6', ...
                                 'dRpH1', 'dRpH2', 'dRpH3', 'dRpH4', 'dRpH5'};
    mono_main('ICE040_IED2') = {'dRpIN2'};
    mono_main('ICE040_IED3') = {'dLaH3', 'dLA5'};

    mono_main('ICE041_IED1') = {'dRTcx4', 'dRTcx5', 'dRTcx6', 'dRTP4', 'dRTP5', 'dRTP6'};
    mono_main('ICE041_IED2') = {'dRmP2', 'dRmP3', 'dRmP4', 'dRmP5', 'dRiP2', 'dRiP3', 'dRiP4', 'dRiP5', 'dRiP6'};
    mono_main('ICE041_IED3') = {'dRPOP2', 'dRPOP3', 'dRPOP4', 'dRPOP5'};
    mono_main('ICE041_IED4') = {'dRlOF3', 'dRlOF4', 'dRlOF5', 'dRlOF6'};

    mono_main('ICE042_IED1') = {'dRpC6', 'dRpC7', 'dRpC8'};
    mono_main('ICE042_IED2') = {'dRmC1', 'dRmC2', 'dRmC3', 'dRSMA1', 'dRSMA2', 'dRSMA3'};
    mono_main('ICE042_IED3') = {'dRmC7', 'dRmC8', 'dRaC5', 'dRaC6'};
    mono_main('ICE042_IED4') = {'dRlOF8', 'dRlOF9', 'dRlOF10'};

    mono_main('ICE043_IED1') = {'dRaH1', 'dRaH2'};
    mono_main('ICE043_IED2') = {'dLaH1', 'dLaH2', 'dLA1', 'dLA2', 'dLA3', 'dLpH1', 'dLpH2'};

    mono_main('ICE044_IED1') = {'dLaH1', 'dLaH2'};
    mono_main('ICE044_IED2') = {'dRaH1', 'dRaH2', 'dRaH3'};
    mono_main('ICE044_IED3') = {'dLpIN3', 'dLpIN4'};

    mono_main('ICE045_IED1') = {'dLaH1', 'dLaH2', 'dLaH3', 'dLpH1', 'dLpH2'};
    mono_main('ICE045_IED2') = {'dRA1', 'dRA2', 'dRA3', 'dRaH1', 'dRaH2', 'dRaH3'};

    mono_main('ICE046_IED1') = {'dRUmTg2'};

    mono_main('ICE047_IED1') = {'dRaH1', 'dRaH2', 'dRaH3'};
    mono_main('ICE047_IED2') = {'dLpH1', 'dLpH2', 'dLpH3'};

    mono_main('ICE048_IED1') = {'dLaH1', 'dLaH2', 'dLpH1', 'dLH2', 'dLpH3'};

    mono_main('ICE049_IED1') = {'dLaH1', 'dLaH2', 'dLaH3', 'dLA1', 'dLA2', 'dLpH1', 'dLpH2', 'dLpH3'};
    mono_main('ICE049_IED2') = {'dLA6', 'dLA7', 'dLaH5', 'dLaH6', 'dLaH7', 'dLasT1', 'dLasT2', 'dLasT3', 'dLasT4', 'dLasT5'};
    mono_main('ICE049_IED3') = {'dLA6', 'dLA7'};

    mono_main('ICE050_IED1') = {'dLM12', 'dLM13', 'dLM14', 'dLSMA3', 'dLcOP3', ...
                                 'dLpIN1', 'dLpIN2', 'dLpIN3', 'dLpIN4', 'dLpIN5', 'dLpIN6'};

    mono_main('ICE051_IED1') = {'dRlmF4', 'dRlmF5', 'dRlmF6', 'dRlmF7', 'dRlmF8', 'dRaIN2', 'dRaIN3'};

    mono_main('ICE052_IED1') = {'dRpsT1', 'dRpsT2', 'dRpsT3', 'dRpsT4'};
    mono_main('ICE052_IED2') = {'dRPOP1', 'dRPOP2', 'dRPOP3', 'dRPOP4'};

    mono_main('ICE053_IED1') = {'dLaH1', 'dLaH2', 'dLaH3', 'dLpH2', 'dLpH3'};

    mono_main('ICE054_IED1') = {'dRiP5', 'dRiP6'};
    mono_main('ICE054_IED2') = {'dRpsT5', 'dRpsT6'};
    mono_main('ICE054_IED3') = {'dRiP5', 'dRiP6', 'dRpsT5', 'dRpsT6'};

    mono_main('ICE055_IED1') = {'dRaH1', 'dRaH2', 'dRaH3', 'dRaH4'};
    mono_main('ICE055_IED2') = {'dLaH1', 'dLaH2'};

    mono_main('ICE056_IED1') = {'dLA3', 'dLaH1', 'dLaH2', 'dLpH1', 'dLpH2'};
    mono_main('ICE056_IED2') = {'dRA1', 'dRA2', 'dRaH1', 'dRaH2', 'dRpH1', 'dRpH2'};

    mono_main('ICE057_IED1') = {'dRaIN1', 'dRaIN2', 'dRaIN3', 'dRaIN4'};
    mono_main('ICE057_IED2') = {'dLmOF5', 'dLmOF6', 'dLmOF7', 'dLmOF8', 'dLmOF9'};

    mono_main('ICE058_IED1') = {'dLaH2', 'dLaH3', 'dLpH1', 'dLpH2', 'dLA1', 'dLA2', 'dLA3'};

    mono_main('ICE059_IED1') = {'dRA5', 'dRA6', 'dRA7', 'dRA8'};
    mono_main('ICE059_IED2') = {'dLA6', 'dLA7', 'dLA8'};
    mono_main('ICE059_IED3') = {'dLA2', 'dLA3'};
    mono_main('ICE059_IED4') = {'dRaH1', 'dRaH2', 'dRA1', 'dRA2'};

    mono_main('ICE060_IED1') = {'dRA1', 'dRA2', 'dRaH1', 'dRaH2', 'dRaH3', 'dRaH4', ...
                                 'dRpH1', 'dRpH2', 'dRpH3', 'dRpH4'};

    mono_main('ICE062_IED1') = {'dRaH1', 'dRaH2','dRA1', 'dRA2', 'dRA3', 'dRpIN4'};
    mono_main('ICE062_IED2') = {'dLaH2', 'dLA1', 'dLA2', 'dLpIN1', 'dLpIN2'};

    mono_main('ICE063_IED1') = {'dLpH1', 'dLpH2', 'dLaH2', 'dLaH3'};
    mono_main('ICE063_IED2') = {'dLA1', 'dLA2'};

    mono_main('ICE064_IED1') = {'dRaH1', 'dRaH2', 'dRaH3', 'dRaH4'};
    mono_main('ICE064_IED2') = {'dRsO5', 'dRsO6', 'dRsO7', 'dRsO8'};

    mono_main('ICE065_IED1') = {'dLpH1', 'dLpH2', 'dLpH3', 'dLaH1', 'dLaH2', 'dLaH3'};

    mono_main('ICE066_IED1') = {'dRpIN1', 'dRpIN2', 'dRpIN3'};
    mono_main('ICE066_IED2') = {'dRaH1', 'dRaH2', 'dRpH2'};

    mono_main('ICE069_IED1') = {'dLaH1', 'dLaH2', 'dLA1', 'dLA2'};
    mono_main('ICE069_IED2') = {'dRaH1', 'dRaH2', 'dRaH3', 'dRA1', 'dRA2'};

    mono_main('ICE070_IED1') = {'dRaIN1', 'dRaIN2'};
    mono_main('ICE070_IED2') = {'dRasT3', 'dRasT4', 'dRasT5'};
    mono_main('ICE070_IED3') = {'dRPOP4', 'dRPOP5'};

    % ---------- 2. Build BIPOLAR main channels using true existing channels ----------
    main_channels_map_bipolar = containers.Map();
    keys_list = mono_main.keys;

    for k = 1:numel(keys_list)
        mapKey = keys_list{k};                 % e.g. 'ICE017_IED1'
        parts  = split(mapKey, '_');
        subjID = parts{1};                     % 'ICE017'

        % --- read existing monopolar channels from ProcessingLog ---
        try
            subjChannels = read_subject_monopolar_channels(baseDataDir, subjID);
        catch ME
            warning('Could not read ProcessingLog for %s: %s. Skipping key %s.', ...
                    subjID, ME.message, mapKey);
            continue;
        end

        subjChannels = unique(subjChannels(:));       % cellstr of monopolar labels
        mainMonoList = mono_main(mapKey);             % cellstr of monopolar main contacts

        % small helper for membership
        isPresent = @(lab) any(strcmp(lab, subjChannels));

        % collect bipolar labels in a small set to avoid duplicates
        bipolarSet = containers.Map('KeyType','char', 'ValueType','logical');

        for m = 1:numel(mainMonoList)
            thisMono = char(mainMonoList{m});   % e.g. 'dRaIN3'
            [prefix, idx0] = split_prefix_num(thisMono);

            if isnan(idx0)
                warning('Cannot parse index from main channel "%s" (key %s). Skipping.', thisMono, mapKey);
                continue;
            end

            % Look at neighbors: idx0-1 and idx0+1
            for delta = [-1, +1]
                nbrIdx = idx0 + delta;
                if nbrIdx <= 0
                    continue;
                end

                nbrMono = sprintf('%s%d', prefix, nbrIdx);

                if isPresent(nbrMono)
                    lo = min(idx0, nbrIdx);
                    hi = max(idx0, nbrIdx);
                    bipLabel = sprintf('%s%d-%s%d', prefix, lo, prefix, hi);
                    bipolarSet(bipLabel) = true;
                end
            end
        end

        bipKeys = bipolarSet.keys;
        if isempty(bipKeys)
            warning('No valid bipolar main channels found for key %s (subject %s).', mapKey, subjID);
            continue;
        end

        % sort for reproducibility
        bipList = sort(bipKeys);

        main_channels_map_bipolar(mapKey) = bipList;
        fprintf('Built %d bipolar main channels for %s\n', numel(bipList), mapKey);
    end

    % ---------- 3. Save result ----------
    save('main_channels_map_bipolar.mat', 'main_channels_map_bipolar');
    fprintf('Saved main_channels_map_bipolar.mat with %d keys.\n', numel(main_channels_map_bipolar.keys));

end

% =====================================================================
function labels = read_subject_monopolar_channels(baseDataDir, subjID)
% Reads the subject's monopolar channel labels from:
% <baseDataDir>/<subjID>/3_EEG/2_Cleaned/<subjID>_Run1_Cleaned/IED_Cleaned/<subjID>_Run1_ProcessingLog_&_ChannelLabels.txt

    logDir  = fullfile(baseDataDir, subjID, '3_EEG', '2_Cleaned', ...
                       [subjID '_Run1_Cleaned'], 'IED_Cleaned');
    logFile = fullfile(logDir, [subjID '_Run1_ProcessingLog_&_ChannelLabels.txt']);

    if ~isfile(logFile)
        error('ProcessingLog file not found: %s', logFile);
    end

    txt = fileread(logFile);

    % Extract block between 'c_Labels (Row):' and 'c_Labels (Column):'
    startIdx = strfind(txt, 'c_Labels (Row):');
    if isempty(startIdx)
        error('Could not find "c_Labels (Row):" in %s', logFile);
    end

    colIdx = strfind(txt, 'c_Labels (Column):');
    if isempty(colIdx)
        segment = txt(startIdx:end);
    else
        segment = txt(startIdx:colIdx-1);
    end

    % Drop the first line containing 'c_Labels (Row):'
    newlinePos = strfind(segment, sprintf('\n'));
    if ~isempty(newlinePos)
        segment = segment(newlinePos(1)+1:end);
    end

    % Replace newlines/tabs with spaces, then split
    segment(segment == sprintf('\n')) = ' ';
    segment(segment == sprintf('\r')) = ' ';
    segment(segment == sprintf('\t')) = ' ';

    tokens = strsplit(strtrim(segment));
    labels = tokens(~cellfun(@isempty, tokens));
end

% =====================================================================
function [prefix, idx] = split_prefix_num(label)
% Splits a label like 'dRaIN3' into:
%   prefix = 'dRaIN'
%   idx    = 3

    label = char(label);
    m = regexp(label, '^([A-Za-z]+)(\d+)$', 'tokens', 'once');

    if isempty(m)
        prefix = label;
        idx    = NaN;
    else
        prefix = m{1};
        idx    = str2double(m{2});
    end
end
