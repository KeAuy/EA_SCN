%RESIDUALGEOMETRIC Finite-difference residual on the geometric grid.
%
%   Discretizes lambda (b'' + 2 b'/r) = (b - 1)(nu b + eta) + b/r with central
%   differences, b(1) = b0 at the inner boundary and the two-term far-field
%   expansion at the outer boundary.

function F = residualGeometric(b, r, hLeft, hRight, b0, lambda, nu, eta)

n = numel(b);
F = zeros(n, 1);

F(1) = b(1) - b0;

for k = 2:n-1
    hMean = 0.5 * (hLeft(k) + hRight(k));

    F(k) = lambda * ((b(k+1) - b(k)) / hRight(k) ...
                   - (b(k) - b(k-1)) / hLeft(k)) / hMean ...
         + 2 * lambda * (b(k+1) - b(k-1)) ...
           / ((hLeft(k) + hRight(k)) * r(k)) ...
         - (b(k) - 1) * (nu * b(k) + eta) ...
         - b(k) / r(k);
end

F(n) = b(n) - buffer.farField(r(n), lambda, nu, eta, 2);
end
