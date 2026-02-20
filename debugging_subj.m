subj='ICE013';
C = load(['/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/' subj '/config_' subj '.mat'],'C'); C=C.C;
EEG = load(fullfile(C.outDir, ['EEG_segments_' subj '.mat']), 'EEG'); EEG=EEG.EEG;

runKey='Run1a';
monoLabels = string(EEG.(runKey).labels);
monoLabels = upper(strrep(monoLabels,'_','-'));

% show anything that looks like DLA
disp(monoLabels(contains(monoLabels,'DLA')))

% check specifically for DLA4 variants
disp(monoLabels(contains(monoLabels,'DLA4') | contains(monoLabels,'DLA-4') | contains(monoLabels,'DLA04') | contains(monoLabels,'DLA 4')))

%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% checking for a specific channel like 'DLA4'
subj='ICE013';
C = load(['/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/' subj '/config_' subj '.mat'],'C'); C=C.C;
EEG = load(fullfile(C.outDir, ['EEG_segments_' subj '.mat']), 'EEG'); EEG=EEG.EEG;

runKey='Run1a';
monoLabels = string(EEG.(runKey).labels);
monoLabels = upper(strrep(monoLabels,'_','-'));

% show anything that looks like DLA
disp(monoLabels(contains(monoLabels,'DLA')))

% check specifically for DLA4 variants
disp(monoLabels(contains(monoLabels,'DLA4') | contains(monoLabels,'DLA-4') | contains(monoLabels,'DLA04') | contains(monoLabels,'DLA 4')))
