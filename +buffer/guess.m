%GUESS Initial guess for BVP5C built from the two series expansions.
%
%   Uses the near-origin branch where the two expansions differ by less than
%   10 percent (or where r < 1) and the far-field branch elsewhere.

function y0 = guess(r, b0, lambda, nu, eta)

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
