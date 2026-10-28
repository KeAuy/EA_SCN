function b = doubleExpSeries(r, lambda, nu, eta, realBranchIndex)
%DOUBLEEXPSERIES Double exponential series of the buffer solution.
%
%   B = BUFFER.DOUBLEEXPSERIES(R, LAMBDA, NU, ETA, REALBRANCHINDEX) uses the
%   cubic coefficients P = 2 lambda/3, Q = lambda - eta q^2/2 and R = 1.
%
%   REALBRANCHINDEX > 0 makes the columns with q < q(REALBRANCHINDEX) real,
%   reproducing the branch switch of Figures 2/4/5.  Index 0 (default) keeps
%   the complex profile for every column, as in Figure 1.
%
%   See also BUFFER.DOUBLEEXPALPHA, BUFFER.DOUBLEEXPPROFILE.

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
