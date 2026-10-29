%JACOBIANGEOMETRICLAMBDA Tridiagonal Jacobian for every lambda column.
%
%   [LOWER, MAIN, UPPER] each have one column per lambda and are laid out for
%   BUFFER.SOLVETRIDIAGONALBATCH: LOWER(1) and UPPER(end) are unused.

function [lower, main, upper] = jacobianGeometricLambda(b, r, hLeft, hRight, lambdaValues, nu, eta)
arguments
    b double
    r double
    hLeft double
    hRight double
    lambdaValues double
    nu (1,1) double
    eta (1,1) double
end

[n, nLambda] = size(b);
lambda       = lambdaValues(:).';

if numel(lambda) ~= nLambda
    error('buffer:sizeMismatch', ...
        'LAMBDA must have one entry per column of B.');
end

lower = zeros(n, nLambda);
upper = zeros(n, nLambda);
main  = zeros(n, nLambda);

k     = (2:n-1).';
hMean = 0.5 * (hLeft(k) + hRight(k));
drSum = hLeft(k) + hRight(k);
rk    = r(k);

for i = 1:nLambda
    lam = lambda(i);

    lower(k, i) = lam ./ hLeft(k) ./ hMean ...
        - 2 * lam ./ (drSum .* rk);
    upper(k, i) = lam ./ hRight(k) ./ hMean ...
        + 2 * lam ./ (drSum .* rk);
    main(2:n-1, i) = (-lam ./ hRight(k) - lam ./ hLeft(k)) ./ hMean ...
        - (2 * nu * b(k, i) - nu + eta) - 1 ./ rk;
end

main(1, :) = 1;
main(n, :) = 1;
end
