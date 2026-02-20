%% written by Tahereh Rashnavadi
% Jan 2026

%%

function [bipData, bipNames] = build_bipolar_signals_from_monopolar(monoData, monoLabels)
% Build bipolar montage from monopolar depth contacts by shaft.
% bipolar = contact(k) - contact(k+1) for adjacent contacts on same shaft.
%
% INPUT:
%   monoData   [nChan x nSamp]
%   monoLabels [nChan x 1] string (e.g., dRHG3, dRHG4)
% OUTPUT:
%   bipData    [nBip x nSamp]
%   bipNames   [nBip x 1] string (e.g., dRHG3-dRHG4)

monoLabels = string(monoLabels);

% Parse prefix + number
prefix = regexprep(monoLabels, '\d+$', '');
numStr = regexp(monoLabels, '\d+$', 'match', 'once');
numVal = nan(size(monoLabels));
for i = 1:numel(monoLabels)
    if ~isempty(numStr{i})
        numVal(i) = str2double(numStr{i});
    end
end

% Keep only labels that have numeric suffix
valid = ~isnan(numVal);
monoData = monoData(valid,:);
monoLabels = monoLabels(valid);
prefix = prefix(valid);
numVal = numVal(valid);

u = unique(prefix);

bipData = [];
bipNames = strings(0,1);

for s = 1:numel(u)
    p = u(s);
    idx = find(prefix == p);
    [~,ord] = sort(numVal(idx));
    idx = idx(ord);
    nums = numVal(idx);

    % adjacent pairs
    for k = 1:numel(idx)-1
        if nums(k+1) == nums(k) + 1
            i1 = idx(k);
            i2 = idx(k+1);

            % bipolar = first - second (consistent convention)
            bip = monoData(i1,:) - monoData(i2,:);

            lab1 = monoLabels(i1);
            lab2 = monoLabels(i2);
            bipLab = lab1 + "-" + lab2;

            bipData = [bipData; bip]; %#ok<AGROW>
            bipNames(end+1,1) = bipLab; %#ok<AGROW>
        end
    end
end
end
