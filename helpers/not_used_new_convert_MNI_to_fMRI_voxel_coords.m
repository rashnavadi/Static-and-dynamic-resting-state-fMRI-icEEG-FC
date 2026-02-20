%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch std2imgcoord from MNI → functional Space for multiple subjects
% 
% This script:
%  1) Loops over a list of subjects.
%  2) For each subject, reads an Excel file containing coordinates in MNI space.
%  3) Extracts x,y,z columns and writes them to a temporary text file.
%  4) Calls `std2imgcoord` to transform MNI coords → native space (example_func).
%  5) Reads back the transformed coordinates, places them into the table,
%     and writes out a new Excel file with updated coordinates.
%
% Written by: Tahereh Rashnavadi, Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Define directories and subject list
ICE_reg_folder       = '/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/';
input_coords_folder  = '/Volumes/MIND/ICE/original_ICE/MNI_coordinates';
output_coords_folder = '/Volumes/MIND/ICE/Tara/native_space_electrodes_coords';

% Make sure the output folder exists
if ~exist(output_coords_folder, 'dir')
    mkdir(output_coords_folder);
end

% Define subjects - for example, ICE001 through ICE070
subject_list = arrayfun(@(n) sprintf('ICE%03d', n), 1:70, 'UniformOutput', false);

%% 2) Loop over subjects
for iSubj = 13:length(subject_list)   
    subject_id = subject_list{iSubj};
    fprintf('\n=== Processing %s ===\n', subject_id);
    
    % Path to subject's Excel file (with MNI coords)
    excel_file = fullfile(input_coords_folder, [subject_id '_channel_info.xlsx']);
    
    % Skip if the file doesn't exist
    if ~isfile(excel_file)
        fprintf('  --> No Excel file found for %s. Skipping.\n', subject_id);
        continue;
    end

    % Read the full Excel table
    T = readtable(excel_file);

    % Check for required columns
    requiredCols = {'x', 'y', 'z', 'ElectrodeType'};
    if ~all(ismember(requiredCols, T.Properties.VariableNames))
        fprintf('Skipping %s: Missing required columns in electrode file.\n', subject_id);
        continue;
    end

    %% Convert x, y, z columns to numeric if they aren't already
    % (e.g., if they're stored as strings in Excel)
    if ~isnumeric(T.x)
        T.x = str2double(T.x);
    end
    if ~isnumeric(T.y)
        T.y = str2double(T.y);
    end
    if ~isnumeric(T.z)
        T.z = str2double(T.z);
    end

    % Handle any non-numeric entries that couldn't be converted
    T.x(isnan(T.x)) = NaN;
    T.y(isnan(T.y)) = NaN;
    T.z(isnan(T.z)) = NaN;

    %% Identify depth electrodes
    is_depth = strcmpi(T.ElectrodeType, 'depth');
    if sum(is_depth) == 0
        fprintf('No depth electrodes found for %s. Skipping.\n', subject_id);
        continue;
    end

    % Extract only the MNI coords for depth electrodes
    mni_coords = T{is_depth, {'x', 'y', 'z'}};

    % If all MNI coords are NaN or empty for depth, skip
    if all(isnan(mni_coords(:)))
        fprintf('All depth coords are NaN for %s. Skipping.\n', subject_id);
        continue;
    end

    % Create a temporary file for writing MNI coords (for FSL)
    mni_coords_file = fullfile(output_coords_folder, sprintf('%s_MNI_coords.txt', subject_id));
    writematrix(mni_coords, mni_coords_file, 'Delimiter', 'space');

    % Define output file for transformed native-space coordinates
    native_coords_file = fullfile(output_coords_folder, sprintf('%s_native_coords.txt', subject_id));

    %----------------------------------------------------------------------
    % 3) Find all runs (Run*, e.g. Run1, Run1a, Run2, Run3b, etc.)
    %----------------------------------------------------------------------
    % --- Correct subject preprocessing folder ---
    preproc_dir = fullfile(ICE_reg_folder, subject_id, ...
        '4_MRI', '2_Functionals', '2_PreProcessing');

    if ~isfolder(preproc_dir)
        fprintf('  --> PreProcessing folder not found for %s: %s\n', subject_id, preproc_dir);
        continue;
    end

    % --- Find Run*.feat folders, excluding "*Exclude*" ---
    run_feats = dir(fullfile(preproc_dir, 'Run*.feat'));
    run_feats = run_feats([run_feats.isdir]);

    % remove Exclude
    keep = ~contains({run_feats.name}, 'Exclude', 'IgnoreCase', true);
    run_feats = run_feats(keep);

    if isempty(run_feats)
        fprintf('  --> No usable Run*.feat (non-Exclude) found for %s in %s\n', subject_id, preproc_dir);
        continue;
    end

    % --- Loop through runs and use reg folder ---
    for iRun = 1:numel(run_feats)
        run_feat_name = run_feats(iRun).name;   % e.g., "Run1a.feat"
        reg_dir = fullfile(preproc_dir, run_feat_name, 'reg');

        fprintf('\n  --> Checking %s / %s\n', subject_id, run_feat_name);

        example_func_nii = fullfile(reg_dir, 'example_func.nii.gz');
        standard_nii     = fullfile(reg_dir, 'standard.nii.gz');
        xfm_mat          = fullfile(reg_dir, 'example_func2standard.mat');

        if ~exist(example_func_nii,'file') || ~exist(standard_nii,'file') || ~exist(xfm_mat,'file')
            fprintf('  --> Missing example_func/standard/mat in %s. Skipping this run.\n', reg_dir);
            continue;
        end

        %----------------------------------------------------------------------
        % 5) Build std2imgcoord command:
        %    standard(=MNI), image(=example_func), and the matrix
        %    that goes from example_func → standard 
        % NOTE:
        % In FSL, voxel coordinates are zero-indexed (i.e., the first voxel has coordinate (0,0,0)),
        % whereas in MATLAB, indices start from 1. This means that if you want to access the transformed
        % voxel coordinates in MATLAB, you need to add 1 to each coordinate obtained from std2imgcoord.
        % MATLAB uses one-based indexing, meaning that voxel (0,0,0) in FSL is actually (1,1,1) in MATLAB.
        % if you extract timeseries in matlab you need to Add 1 to shift
        % from FSL's zero-based index to MATLAB's one-based index,
        %----------------------------------------------------------------------
        tmp_native_txt = fullfile(output_coords_folder, [subject_id '_native_tmp_mm.txt']);

        % Output text file (transformed coords from MNI space to func
        % space), Output will be 0-based voxel coordinates in example_func space.
        fsl_cmd = sprintf([ ...
            'std2imgcoord ' ... % example_func image (target space)
            '-img "%s" ' ...    % standard/MNI  (source space)
            '-std "%s" ' ...    % example_func->standard matrix
            '-xfm "%s" ' ...    % input coords -> output coords
            '-mm ' ...          % <-- REQUIRED (input is MNI mm)
            '%s > %s'], ...
            example_func_nii, ...
            standard_nii, ...
            xfm_mat, ...
            mni_coords_file, ...
            tmp_native_txt);

        fprintf('  --> Running FSL: %s\n', fsl_cmd);

        [status, cmdout] = system(fsl_cmd);
        if status ~= 0
            fprintf('  --> std2imgcoord failed for %s:\n%s\n', subject_id, cmdout);
            continue;
        end

        %----------------------------------------------------------------------
        % 6) Read back the transformed coords and put them in the table
        %----------------------------------------------------------------------
        if ~isfile(tmp_native_txt)
            fprintf('  --> No transformed file output for %s. Skipping.\n', subject_id);
            continue;
        end

        native_coords_mm = readmatrix(tmp_native_txt);  % mm in example_func space

        if size(native_coords_mm, 2) ~= 3
            fprintf('  --> Transformed coords do not have 3 columns for %s. Skipping.\n', subject_id);
            continue;
        end

        % ---------------------------------------------------------
        % Convert mm -> voxel using example_func affine
        % ---------------------------------------------------------
        V = niftiread(example_func_nii); %#ok<NASGU>
        info = niftiinfo(example_func_nii);
        A = info.Transform.T';   % affine matrix (voxel -> mm)

        % Convert mm -> voxel:  ijk = inv(A) * [x y z 1]'
        Ainv = inv(A);

        nPts = size(native_coords_mm,1);
        xyz1 = [native_coords_mm, ones(nPts,1)];
        vox0 = (Ainv * xyz1')';
        vox0 = vox0(:,1:3);                % continuous 0-based vox
        vox0_round = round(vox0);          % integer 0-based vox

        fprintf("Voxel output ranges: min=[%.2f %.2f %.2f], max=[%.2f %.2f %.2f]\n", ...
            min(vox0_round,[],1), max(vox0_round,[],1));

        % ---------------------------------------------------------
        % Save voxel coordinates into table
        % ---------------------------------------------------------
        T{is_depth, {'x','y','z'}} = vox0_round;
        
        %----------------------------------------------------------------------
        % 7) Save a new Excel file with the updated (native) coords
        %----------------------------------------------------------------------
        out_excel = fullfile(output_coords_folder, [subject_id '_fmri_space_vox_coords.xlsx']);
        writetable(T, out_excel, 'FileType', 'spreadsheet');
        fprintf('  --> Saved transformed coords for %s: %s\n', subject_id, out_excel);
        
        break;   % stop after first successful run and stop to loop through all runs

    end
    % *Always* delete the temporary files at the end of the loop
    if exist(mni_coords_file, 'file'), delete(mni_coords_file); end
    if exist(tmp_native_txt, 'file'), delete(tmp_native_txt); end
end


fprintf('\nAll done!\n');
