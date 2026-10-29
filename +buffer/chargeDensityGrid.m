function c = chargeDensityGrid(b, r, nuValues)
%CHARGEDENSITYGRID Charge density on a lambda-nu grid.
%
%   C = BUFFER.CHARGEDENSITYGRID(B, R, NUVALUES) evaluates
%
%       c = nu (b - 1) + 1/r
%
%   for B of size NR-by-NLAMBDA-by-NNU.  NUVALUES is the vector of nu values
%   along the third dimension.
%
%   See also BUFFER.CHARGEDENSITY, BUFFER.MEANERRORGRID.

r  = reshape(r, [], 1);
nu = reshape(nuValues, 1, 1, []);

if size(b, 3) ~= numel(nuValues)
    error('buffer:sizeMismatch', ...
        'The third dimension of B must match NUMEL(NUVALUES).');
end

c = nu .* (b - 1) + 1 ./ r;
end
