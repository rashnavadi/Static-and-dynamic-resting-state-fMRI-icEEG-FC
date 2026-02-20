%% written by Tahereh Rashnavadi
% Jan 2026build_bipolar_signals_from_monopolar
function out = canonical_bipolar_name(in)
%CANONICAL_BIPOLAR_NAME Standardize bipolar labels to "PREFIXlo-PREFIXhi".
%
% Examples:
%   "dLA4-dLA3"   -> "DLA3-DLA4"
%   "dLA4_ dLA3"  -> "DLA3-DLA4"
%   "dRaIN3-dRaIN4" -> "DRAIN3-DRAIN4"
%   "A1-B2"       -> "A1-B2"  (prefix differs; preserved order, just uppercase)
%
% Rules:
% - Converts '_' to '-' and removes whitespace
% - Uppercases output
% - If both sides share the same prefix and both end with integers, sort indices ascending
% - Otherwise, returns the cleaned uppercase input unchanged (but with '-' delimiter)

    in = string(in);

    % Normalize delimiters and whitespace
    s = upper(in);
    s = replace(s, " ", "");
    s = replace(s, "_", "-");

    out = strings(size(s));

    for i = 1:numel(s)
        si = s(i);

        % Split on dash
        parts = split(si, "-");
        if numel(parts) ~= 2
            out(i) = si;  % not a simple bipolar label
            continue;
        end

        a = parts(1);
        b = parts(2);

        [pa, ia, okA] = split_prefix_num(a);
        [pb, ib, okB] = split_prefix_num(b);

        % If we can't parse, just return cleaned label
        if ~okA || ~okB
            out(i) = si;
            continue;
        end

        % If prefixes differ, keep order (but standardized)
        if ~strcmp(pa, pb)
            out(i) = pa + string(ia) + "-" + pb + string(ib);
            continue;
        end

        % Same prefix: sort indices
        lo = min(ia, ib);
        hi = max(ia, ib);
        out(i) = pa + string(lo) + "-" + pa + string(hi);
    end
end

function [prefix, idx, ok] = split_prefix_num(label)
% Split something like "DRAIN3" into prefix="DRAIN", idx=3
% Returns ok=false if it doesn't match.

    label = char(label);
    m = regexp(label, '^([A-Z]+)(\d+)$', 'tokens', 'once');
    if isempty(m)
        prefix = string(label);
        idx = NaN;
        ok = false;
    else
        prefix = string(m{1});
        idx = str2double(m{2});
        ok = isfinite(idx);
    end
end
