%CHARGEDENSITY Charge-density profile c = nu (b - 1) + 1/r.
%
%   B is NUMEL(R)-by-NUMEL(NU); the result has the same size.

function c = chargeDensity(b, r, nu)

r  = r(:);
nu = nu(:)';
c  = nu .* (b - 1) + 1 ./ r;
end
