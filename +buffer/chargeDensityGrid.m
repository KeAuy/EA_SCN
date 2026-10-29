%CHARGEDENSITYGRID Charge density on a lambda-nu grid.
%
%   Evaluates c = nu (b - 1) + 1/r for B of size NR-by-NLAMBDA-by-NNU, with
%   NUVALUES running along the third dimension.

function c = chargeDensityGrid(b, r, nuValues)

r  = reshape(r, [], 1);
nu = reshape(nuValues, 1, 1, []);

if size(b, 3) ~= numel(nuValues)
    error('buffer:sizeMismatch', ...
        'The third dimension of B must match NUMEL(NUVALUES).');
end

c = nu .* (b - 1) + 1 ./ r;
end
