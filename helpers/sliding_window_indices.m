%% written by Tahereh Rashnavadi
% Jan 2026build_bipolar_signals_from_monopolar

function [winStarts, winEnds] = sliding_window_indices(T, winLen, stepLen, lag)
%SLIDING_WINDOW_INDICES Generate sliding window indices with optional lag shift.
%
% Inputs:
%   T       = total length in samples/volumes
%   winLen  = window length (samples/volumes)
%   stepLen = step size (samples/volumes)
%   lag     = shift forward in samples/volumes (applied to windows)
%
% Outputs:
%   winStarts, winEnds = vectors of window start/end indices (1-based, inclusive)

    if nargin < 4, lag = 0; end
    if any([T, winLen, stepLen] <= 0)
        error('Inputs T, winLen, stepLen must be positive.');
    end
    if winLen > T
        winStarts = [];
        winEnds = [];
        return;
    end

    % base windows (no lag) in 1..T
    baseStarts = 1:stepLen:(T - winLen + 1);
    baseEnds   = baseStarts + winLen - 1;

    % apply lag shift (forward)
    winStarts = baseStarts + lag;
    winEnds   = baseEnds   + lag;

    % keep only windows fully inside 1..T
    ok = (winStarts >= 1) & (winEnds <= T);
    winStarts = winStarts(ok);
    winEnds   = winEnds(ok);
end
