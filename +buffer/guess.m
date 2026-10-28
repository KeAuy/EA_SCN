function y0 = guess(r, b0, lambda, nu, eta)
%GUESS Initial guess for BVP5C built from the two series expansions.
%
%   Y0 = BUFFER.GUESS(R, B0, LAMBDA, NU, ETA) returns the state [b; db/dr]
%   at the radii R.  The near-origin branch is used where the two expansions
%   agree to within 10% and for r < 1; the far-field branch is used
%   elsewhere.
%
%   See also BUFFER.EXPANSION, BUFFER.ODE, BUFFER.BC.

r = r(:);
[bNear, dbNear, bFar, dbFar] = buffer.expansion(r, b0, lambda, nu, eta, 4);

agreement = abs(bNear - bFar) ./ (abs(bNear) + abs(bFar) + 1e-14);
useNear   = agreement < 0.1 | r < 1;

b  = bFar;
db = dbFar;
b(useNear)  = bNear(useNear);
db(useNear) = dbNear(useNear);

y0 = [b(:).'; db(:).'];
end
