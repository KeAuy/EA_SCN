%DOUBLEEXPSERIES Double-exponential series with constant coefficients.
%
%   Uses P = 2 lambda/3, R = 1 and Q = lambda - eta q^2/2.  When
%   REALBRANCHINDEX is positive, the columns with q below q(REALBRANCHINDEX)
%   are forced onto the real branch, matching the original Figure 2/5 panels.

function b = doubleExpSeries(r, lambda, nu, eta, realBranchIndex)

arguments
    r double {mustBePositive}
    lambda (1,1) double {mustBePositive}
    nu double
    eta (1,1) double
    realBranchIndex (1,1) double {mustBeNonnegative} = 0
end

r  = r(:);
nu = nu(:)';
q  = 1 ./ (eta + nu);

pCoef = 2 * lambda / 3 * ones(size(q));
qCoef = lambda - eta * q.^2 / 2;
rCoef = ones(size(q));

alpha = buffer.doubleExpAlpha(q, eta, pCoef, qCoef, rCoef).';
b     = buffer.doubleExpProfile(r, q, alpha, eta);

if realBranchIndex > 0
    idx  = min(max(round(realBranchIndex), 1), numel(q));
    mask = q < q(idx);
    b(:, mask) = real(b(:, mask));
end
end
