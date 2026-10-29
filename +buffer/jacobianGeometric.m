%JACOBIANGEOMETRIC Exact sparse tridiagonal Jacobian of the FD residual.
%
%   The two Dirichlet rows have a unit diagonal.

function J = jacobianGeometric(b, r, hLeft, hRight, lambda, nu, eta)

n = numel(b);
k = (2:n-1)';

hMean = 0.5 * (hLeft(k) + hRight(k));
drSum = hLeft(k) + hRight(k);
rk    = r(k);

lower = lambda ./ hLeft(k) ./ hMean - 2 * lambda ./ (drSum .* rk);
main  = (-lambda ./ hRight(k) - lambda ./ hLeft(k)) ./ hMean ...
        - (2 * nu * b(k) - nu + eta) - 1 ./ rk;
upper = lambda ./ hRight(k) ./ hMean + 2 * lambda ./ (drSum .* rk);

J = spdiags([[lower; 0; 0], ...
             [1; main; 1], ...
             [0; 0; upper]], -1:1, n, n);
end
