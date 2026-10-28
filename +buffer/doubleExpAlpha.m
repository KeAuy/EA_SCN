function alpha = doubleExpAlpha(q, eta, pCoef, qCoef, rCoef)
%DOUBLEEXPALPHA Decaying root of the double-exponential cubic.
%
%   ALPHA = BUFFER.DOUBLEEXPALPHA(Q, ETA, PCOEF, QCOEF, RCOEF) solves
%
%       eta q^2 P alpha^3 - Q alpha^2 - R alpha + 1/q = 0
%
%   with the three cube-root branches and returns the first branch.  All
%   inputs are column vectors of the same length.
%
%   See also BUFFER.DOUBLEEXPSERIES, BUFFER.DOUBLEEXPVARIANT.

q     = q(:);
pCoef = pCoef(:);
qCoef = qCoef(:);
rCoef = rCoef(:);

w = 3 * eta * q.^2 .* pCoef;
e = qCoef.^2 + w .* rCoef;
f = 1.5 * (qCoef .* rCoef - w ./ q) .* w + qCoef.^3;

cubeRoots = [0.5 + sqrt(3)/2 * 1i, 0.5 - sqrt(3)/2 * 1i, -1];

g = zeros(numel(q), 3);
for m = 1:3
    g(:, m) = cubeRoots(m) ...
        .* (sqrt((f - e.^(3/2)) .* (f + e.^(3/2))) + f).^(1/3);
end

alpha = (1 ./ w) .* (qCoef - g - e ./ g);
alpha = alpha(:, 1);
end
