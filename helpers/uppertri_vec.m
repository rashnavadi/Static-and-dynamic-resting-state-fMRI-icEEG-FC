%% written by Tahereh Rashnavadi
% Jan 2026build_bipolar_signals_from_monopolar

function v = uppertri_vec(M)
%UPPERTRI_VEC Vectorize upper triangle (excluding diagonal) of a square matrix.
%   v = uppertri_vec(M) returns a column vector of M(i,j) for i<j.

    if size(M,1) ~= size(M,2)
        error('uppertri_vec:InputNotSquare', 'Input must be a square matrix.');
    end
    mask = triu(true(size(M)), 1);   % strictly upper triangle
    v = M(mask);
end
