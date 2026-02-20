%% step01_make_segment_table.m
% Written by Tahereh Rashnavadi
% Jan 2026

% DO NOT clear here (breaks the outer loop)
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

% Expect C to exist (set by step00_config.m)
if ~exist("C","var") || ~isstruct(C) || ~isfield(C,"subjectID")
    error("C must be set by the calling script (run step00_config.m first).");
end

addpath(C.codeDir);
addpath(fullfile(C.codeDir, "helpers"));

S = load(C.segmentsMat);   % loads "segments"
segments = S.segments;

assert(isfield(segments, C.subjectID), "segments has no field '%s'", C.subjectID);
segStruct = segments.(C.subjectID);

rows = [];
runFields = string(fieldnames(segStruct));

for k = 1:numel(runFields)
    runID = runFields(k);
    trRange = segStruct.(char(runID));  % [startTR endTR], 0-based inclusive

    startTR = trRange(1);
    endTR   = trRange(2);

    nTR = (endTR - startTR + 1);
    startSec = startTR * C.TR;
    endSec   = (endTR + 1) * C.TR;  % end-exclusive
    durSec   = endSec - startSec;

    rows = [rows; table(C.subjectID, runID, startTR, endTR, nTR, startSec, endSec, durSec)];
end

segTable = rows;
segTable.Properties.VariableNames = ["Subject","RunID","startTR","endTR","nTR","startSec","endSec","durSec"];

disp(segTable);

outFile = fullfile(C.outDir, sprintf("segTable_%s.mat", C.subjectID));
save(outFile, "segTable");
fprintf("[OK] Saved %s\n", outFile);
