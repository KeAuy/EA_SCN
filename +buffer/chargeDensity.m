function c = chargeDensity(b, r, nu)
%CHARGEDENSITY Charge-density profile c = nu (b - 1) + 1/r.
%
%   C = BUFFER.CHARGEDENSITY(B, R, NU) evaluates the quantity used by
%   Figure 5.  B is NUMEL(R)-by-NUMEL(NU), R is a column vector and NU is a
%   vector.
%
%   See also BUFFER.MEANERROR.

r  = r(:);
nu = nu(:)';
c  = nu .* (b - 1) + 1 ./ r;
end
