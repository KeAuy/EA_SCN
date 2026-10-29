%SELECTWINNER Select one slice per grid point according to a method index.
%
%   VALUES is NLAMBDA-by-NNU-by-NMETHOD and WINNER is NLAMBDA-by-NNU; the
%   result is the NLAMBDA-by-NNU slice chosen at every point.

function selected = selectWinner(values, winner)

[nLambda, nNu, nMethod] = size(values);
if ~isequal(size(winner), [nLambda, nNu])
    error('buffer:sizeMismatch', ...
        'WINNER must be NLAMBDA-by-NNU.');
end
if any(winner(:) < 1 | winner(:) > nMethod | winner(:) ~= round(winner(:)))
    error('buffer:invalidIndex', ...
        'WINNER must contain integer method indices.');
end

[iLambda, iNu] = ndgrid(1:nLambda, 1:nNu);
selected = values(sub2ind([nLambda, nNu, nMethod], iLambda, iNu, winner));
end
