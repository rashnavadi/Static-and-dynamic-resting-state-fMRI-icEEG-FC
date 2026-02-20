%% written by Tahereh Rashnavadi
% Jan 2026
%  Space: VOXEL indices

% Call it like:
% extract_bipolar_fMRI_seedTS_excludeIED( ...
%     'ICE017','Run1b', ...
%     '/Volumes/MIND/ICE/ICE_denoised_filtered_funcs', ...
%     '/Volumes/MIND/ICE/Tara/native_space_electrodes_coords', ...
%     '/Volumes/MIND/ICE/scripts/Kristina_paper/helpers/main_channels_map_bipolar.mat', ...
%     '/Volumes/MIND/ICE/Tara/Kristina_paper');


%%
function extract_bipolar_fMRI_seedTS_excludeIED(subjectID, runID, baseFMriDir, baseElecDir, mainMapMatPath, projectRoot)
% extract_bipolar_fMRI_seedTS_excludeIED
%
% Builds fMRI bipolar seeds FIRST, excludes epileptogenic/IED contacts BEFORE ROI extraction,
% then extracts mean BOLD time series for each surviving bipolar pair using 3x3x3 voxel ROIs.
%
% INPUTS
%   subjectID      e.g., 'ICE033'
%   runID          e.g., 'Run2a'
%   baseOutDir     e.g., '/Volumes/MIND/ICE/Tara'  (where fMRI_timeseries/ will be saved)
%   baseFMriDir    e.g., '/Volumes/MIND/ICE/ICE_denoised_filtered_funcs'
%   baseElecDir    e.g., '/Volumes/MIND/ICE/Tara/native_space_electrodes_coords'
%   mainMapMatPath e.g., '/Volumes/MIND/ICE/scripts/IED_marking/main_channels_map_bipolar.mat'
%
% OUTPUTS (saved)
%   <baseOutDir>/fMRI_timeseries/<subjectID>/<runID>/
%       <subjectID>_<runID>_seedTS_bipolar_exclIED.mat
%       <subjectID>_<runID>_seedTS_bipolar_exclIED.csv
%       <subjectID>_<runID>_timeseries_extraction_log.csv
%
% Notes
% - Assumes electrode coords are VOXEL INDICES in native functional space (1-based in Excel/MATLAB),
%   converted to FSL 0-based indexing by subtracting 1.
% - Exclusion is conservative: union of all IED types for that subject.

    
    % --- ensure helpers are on path ---
    codeDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(codeDir, "helpers"));

    %% ---------------------------
    % 0) Basic checks
    %% ---------------------------
    if system('which fsl') ~= 0
        error('FSL not found in PATH. Source FSL before running.');
    end

    fmriFile = fullfile(baseFMriDir, subjectID, runID, sprintf('%s_denoised_filtered_func_%s.nii.gz', subjectID, runID));
    if ~isfile(fmriFile)
        error('fMRI file not found: %s', fmriFile);
    end

    electrodeFile = fullfile(baseElecDir, sprintf('%s_native_space_coords.xlsx', subjectID));
    if ~isfile(electrodeFile)
        error('Electrode coord file not found: %s', electrodeFile);
    end

    if ~isfile(mainMapMatPath)
        error('main_channels_map_bipolar.mat not found: %s', mainMapMatPath);
    end

    if nargin < 6 || isempty(projectRoot)
        projectRoot = '/Volumes/MIND/ICE/Tara/Kristina_paper';
    end

    outDir = fullfile(projectRoot, 'fMRI_timeseries', subjectID, runID);
    if ~exist(outDir, 'dir'), mkdir(outDir); end


    %% ---------------------------
    % 1) Read fMRI dims / TR / voxel size (FSL)
    %% ---------------------------
    [status, infoOutput] = system(sprintf('fslinfo %s', fmriFile));
    if status ~= 0, error('fslinfo failed for %s', fmriFile); end

    infoLines = strsplit(infoOutput, '\n');
    nVolumes = NaN; TR = NaN; voxelSize = nan(1,3); dimXYZ = nan(1,3);

    for ll = 1:numel(infoLines)
        line = strtrim(infoLines{ll});
        toks = regexp(line, '\s+', 'split');
        if numel(toks) < 2, continue; end

        switch toks{1}
            case 'dim1',   dimXYZ(1) = str2double(toks{2});
            case 'dim2',   dimXYZ(2) = str2double(toks{2});
            case 'dim3',   dimXYZ(3) = str2double(toks{2});
            case 'dim4',   nVolumes  = str2double(toks{2});
            case 'pixdim1', voxelSize(1) = str2double(toks{2});
            case 'pixdim2', voxelSize(2) = str2double(toks{2});
            case 'pixdim3', voxelSize(3) = str2double(toks{2});
            case 'pixdim4', TR = str2double(toks{2});
        end
    end

    if any(isnan([dimXYZ nVolumes TR]))
        error('Could not parse fslinfo fully for %s', fmriFile);
    end

    fprintf('fMRI: %s | dims=[%d %d %d] | nVol=%d | TR=%.3f\n', runID, dimXYZ(1),dimXYZ(2),dimXYZ(3), nVolumes, TR);

    inBoundsStart = @(x,y,z) (x>=0 && y>=0 && z>=0 && ...
        (x+2) < dimXYZ(1) && (y+2) < dimXYZ(2) && (z+2) < dimXYZ(3));

    %% ---------------------------
    % 2) Load electrode coords + keep depth only + drop missing coords
    %% ---------------------------
    T = readtable(electrodeFile);

    requiredCols = {'ElectrodeName_labConvention_', 'ContactNumber', 'x', 'y', 'z', 'ElectrodeType'};
    if ~all(ismember(requiredCols, T.Properties.VariableNames))
        error('Expected columns not found in %s', electrodeFile);
    end

    chNames = string(T.ElectrodeName_labConvention_) + string(T.ContactNumber); % e.g. dLA + 1 -> dLA1
    coords  = round([T.x, T.y, T.z]);  % voxel indices (1-based)
    types   = string(T.ElectrodeType);

    valid = all(~isnan(coords),2) & strcmpi(types,'depth');
    chNames = chNames(valid);
    coords  = coords(valid,:);

    if isempty(chNames)
        error('No valid depth electrodes with coords for %s', subjectID);
    end

    %% ---------------------------
    % 3) Build conservative excluded CONTACT set from main_channels_map_bipolar (union across all IED types)
    %% ---------------------------
    S = load(mainMapMatPath);
    if ~isfield(S,'main_channels_map_bipolar')
        error('main_channels_map_bipolar not found inside: %s', mainMapMatPath);
    end
    mainMap = S.main_channels_map_bipolar;

    % keys look like ICE033_IED1, ICE033_IED2, ...
    keysAll = string(mainMap.keys);
    subjKeys = keysAll(startsWith(keysAll, subjectID + "_"));
    fprintf('Found %d IED-type keys for %s in main map.\n', numel(subjKeys), subjectID);

    excludedMono = strings(0,1);
    excludedBip  = strings(0,1);

    for kk = 1:numel(subjKeys)
        bipList = string(mainMap(char(subjKeys(kk))));
        bipList = canonical_bipolar_name(upper(strrep(bipList,'_','-')));
        excludedBip = [excludedBip; bipList(:)];

        % derive monopolar contacts from bipolar labels
        for bb = 1:numel(bipList)
            parts = split(bipList(bb), "-");
            if numel(parts)==2
                excludedMono = [excludedMono; parts(1); parts(2)];
            end
        end
    end

    excludedMono = unique(excludedMono);
    excludedBip  = unique(excludedBip);

    fprintf('Conservative exclusion set for %s:\n', subjectID);
    fprintf('  excluded monopolar contacts: %d\n', numel(excludedMono));
    fprintf('  excluded bipolar labels:     %d\n', numel(excludedBip));

    %% ---------------------------
    % 4) Build candidate bipolar pairs (adjacent contacts on same shaft), exclude if either contact is excluded
    %% ---------------------------
    shaft = regexprep(chNames, '\d+$', '');        % dLA
    cnum  = regexp(chNames, '\d+$', 'match', 'once');
    cnum  = cellfun(@str2double, cellstr(cnum));

    uShaft = unique(shaft);

    seedNames = strings(0,1);
    seedPairs = zeros(0,2);   % indices into chNames/coords (row indices)
    exclReason = strings(0,1);

    for s = 1:numel(uShaft)
        this = uShaft(s);
        idx = find(shaft == this);

        % sort by contact number
        [~, ord] = sort(cnum(idx));
        idx = idx(ord);

        for k = 1:numel(idx)-1
            i1 = idx(k);
            i2 = idx(k+1);

            % label like dLA3-dLA4 (dash style to match iEEG bipolar)
            lbl = sprintf('%s%d-%s%d', this, cnum(i1), this, cnum(i2));
            lbl = canonical_bipolar_name(upper(lbl));

            mono1 = upper(chNames(i1));
            mono2 = upper(chNames(i2));

            % Exclude if EITHER contact is excluded (your requirement)
            if any(excludedMono == mono1) || any(excludedMono == mono2)
                seedNames(end+1,1) = lbl; %#ok<AGROW>
                seedPairs(end+1,:) = [i1 i2]; %#ok<AGROW>
                exclReason(end+1,1) = "excluded_contact"; %#ok<AGROW>
                continue
            end

            % Also exclude if label itself is in excluded bipolar (extra safety)
            if any(excludedBip == lbl)
                seedNames(end+1,1) = lbl; %#ok<AGROW>
                seedPairs(end+1,:) = [i1 i2]; %#ok<AGROW>
                exclReason(end+1,1) = "excluded_bipolar_label"; %#ok<AGROW>
                continue
            end

            seedNames(end+1,1) = lbl; %#ok<AGROW>
            seedPairs(end+1,:) = [i1 i2]; %#ok<AGROW>
            exclReason(end+1,1) = ""; %#ok<AGROW>
        end
    end

    if isempty(seedNames)
        error('No bipolar candidates created for %s %s', subjectID, runID);
    end

    nCandidates = numel(seedNames);

    %% ---------------------------
    % 5) Extract fMRI TS for surviving bipolar pairs only (3x3x3 around each contact, union mask)
    %% ---------------------------
    seedTS = nan(nVolumes, numel(seedNames));
    logRows = {};

    for p = 1:numel(seedNames)
        i1 = seedPairs(p,1);
        i2 = seedPairs(p,2);
        lbl = seedNames(p);

        % If excluded, log and skip ROI extraction entirely
        if strlength(exclReason(p)) > 0
            logRows(end+1,:) = {char(lbl), char(chNames(i1)), char(chNames(i2)), ...
                'SKIP_EXCLUDED', char(exclReason(p)), NaN}; %#ok<AGROW>
            continue
        end

        % convert 1-based voxel coords -> 0-based for FSL ROI start
        x1 = round(coords(i1,1)) - 1;  y1 = round(coords(i1,2)) - 1;  z1 = round(coords(i1,3)) - 1;
        x2 = round(coords(i2,1)) - 1;  y2 = round(coords(i2,2)) - 1;  z2 = round(coords(i2,3)) - 1;

        if ~inBoundsStart(x1,y1,z1) || ~inBoundsStart(x2,y2,z2)
            logRows(end+1,:) = {char(lbl), char(chNames(i1)), char(chNames(i2)), ...
                'SKIP_OOB', 'roi_start_out_of_bounds', 0}; %#ok<AGROW>
            continue
        end

        mask1 = fullfile(outDir, sprintf('roi_%s.nii.gz', char(chNames(i1))));
        mask2 = fullfile(outDir, sprintf('roi_%s.nii.gz', char(chNames(i2))));
        pairMask = fullfile(outDir, sprintf('roi_%s.nii.gz', char(strrep(lbl,'-','_'))));
        tsFile = fullfile(outDir, sprintf('%s_%s_%s_bipolar_ts.txt', subjectID, runID, char(strrep(lbl,'-','_'))));

        % build 3x3x3 masks for both contacts
        [sA,oA] = system(sprintf('fslmaths %s -roi %d 3 %d 3 %d 3 0 -1 %s', fmriFile, x1,y1,z1, mask1));
        [sB,oB] = system(sprintf('fslmaths %s -roi %d 3 %d 3 %d 3 0 -1 %s', fmriFile, x2,y2,z2, mask2));
        if sA~=0 || sB~=0
            logRows(end+1,:) = {char(lbl), char(chNames(i1)), char(chNames(i2)), ...
                'FAIL_ROI', strtrim([oA ' | ' oB]), 0}; %#ok<AGROW>
            safeDelete(mask1); safeDelete(mask2); safeDelete(pairMask);
            continue
        end

        % union mask
        system(sprintf('fslmaths %s -add %s %s', mask1, mask2, pairMask));

        % check voxels
        [sV,oV] = system(sprintf('fslstats %s -V', pairMask));
        voxCount = 0; if sV==0, voxCount = str2double(strtok(oV)); end
        if voxCount<=0
            logRows(end+1,:) = {char(lbl), char(chNames(i1)), char(chNames(i2)), ...
                'FAIL_EMPTY', '0 voxels in union mask', voxCount}; %#ok<AGROW>
            safeDelete(mask1); safeDelete(mask2); safeDelete(pairMask);
            continue
        end

        % extract mean TS
        [sM, oM] = system(sprintf('fslmeants -i %s -m %s -o %s', fmriFile, pairMask, tsFile));
        if sM~=0
            logRows(end+1,:) = {char(lbl), char(chNames(i1)), char(chNames(i2)), ...
                'FAIL_MEANTS', strtrim(oM), voxCount}; %#ok<AGROW>
            safeDelete(mask1); safeDelete(mask2); safeDelete(pairMask);
            continue
        end

        ts = load(tsFile);
        if numel(ts) ~= nVolumes || all(~isfinite(ts)) || std(ts(~isnan(ts)))==0
            logRows(end+1,:) = {char(lbl), char(chNames(i1)), char(chNames(i2)), ...
                'FAIL_TS', sprintf('len=%d expected=%d', numel(ts), nVolumes), voxCount}; %#ok<AGROW>
            safeDelete(mask1); safeDelete(mask2); safeDelete(pairMask);
            continue
        end

        seedTS(:,p) = ts(:);
        logRows(end+1,:) = {char(lbl), char(chNames(i1)), char(chNames(i2)), ...
            'OK', '', voxCount}; %#ok<AGROW>

        safeDelete(mask1); safeDelete(mask2); safeDelete(pairMask);
    end

    % Keep only successfully extracted seeds
    ok = all(isfinite(seedTS),1);
    seedTS_bipolar = seedTS(:, ok);
    seedNames_bipolar = seedNames(ok);

    % compute counts from the log to separate IED-excluded vs technical skips
    logT = cell2table(logRows, 'VariableNames', {'BipolarLabel','Contact1','Contact2','Status','Detail','Voxels'});

    nOK        = sum(strcmp(logT.Status, "OK"));
    nSkipIED   = sum(startsWith(logT.Status, "SKIP_EXCLUDED"));  % excluded for IED rule
    nSkipTech  = sum(startsWith(logT.Status, "SKIP_")) - nSkipIED; % e.g., SKIP_OOB
    nFailTech  = sum(startsWith(logT.Status, "FAIL_"));

    fprintf('SUMMARY %s %s:\n', subjectID, runID);
    fprintf('  candidates (adjacent bipolars): %d\n', nCandidates);
    fprintf('  excluded due to IED contacts:   %d\n', nSkipIED);
    fprintf('  skipped technical (e.g., OOB):  %d\n', nSkipTech);
    fprintf('  failed extraction technical:    %d\n', nFailTech);
    fprintf('  final usable bipolar seeds OK:  %d\n', nOK);

    fprintf('DONE %s %s: extracted %d bipolar seeds (excluded/skipped %d).\n', ...
        subjectID, runID, numel(seedNames_bipolar), sum(~ok));

    %% ---------------------------
    % 6) Save outputs
    %% ---------------------------
    outMat = fullfile(outDir, sprintf('%s_%s_seedTS_bipolar_exclIED.mat', subjectID, runID));
    outCsv = fullfile(outDir, sprintf('%s_%s_seedTS_bipolar_exclIED.csv', subjectID, runID));
    outLog = fullfile(outDir, sprintf('%s_%s_timeseries_extraction_log.csv', subjectID, runID));

    save(outMat, 'seedTS_bipolar', 'seedNames_bipolar', 'subjectID', 'runID', 'TR', 'nVolumes', 'voxelSize', ...
    'nCandidates','nOK','nSkipIED','nSkipTech','nFailTech', '-v7.3');


    if ~isempty(seedTS_bipolar)
        writetable(array2table(seedTS_bipolar, 'VariableNames', cellstr(seedNames_bipolar)), outCsv);
    else
        warning('No valid seeds extracted -> CSV not written.');
    end

    if ~isempty(logRows)
        logT = cell2table(logRows, 'VariableNames', {'BipolarLabel','Contact1','Contact2','Status','Detail','Voxels'});
        writetable(logT, outLog);
    end

    summaryT = table(string(subjectID), string(runID), nCandidates, nSkipIED, nSkipTech, nFailTech, nOK, ...
    'VariableNames', {'Subject','Run','nCandidates','nSkipIED','nSkipTech','nFailTech','nOK'});

    writetable(summaryT, fullfile(outDir, sprintf('%s_%s_seed_exclusion_summary.csv', subjectID, runID)));

end


function [prefix, idx] = split_prefix_num(label)
    label = char(label);
    m = regexp(label, '^([A-Za-z]+)(\d+)$', 'tokens', 'once');
    if isempty(m)
        prefix = label; idx = NaN;
    else
        prefix = m{1}; idx = str2double(m{2});
    end
end

function safeDelete(f)
    if exist(f,'file'), delete(f); end
end