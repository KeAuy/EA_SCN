%RESIDUALJACOBIANUNIFORM Residual and Jacobian on an equidistant grid.
%
%   Used by BUFFER.SOLVEUNIFORM for the Figure 1 and 3 panels.  The outer row
%   uses the four-term far-field expansion.

function [F, J] = residualJacobianUniform(b, r, b0, lambda, nu, eta)

n = numel(b);
h = r(2) - r(1);

F = zeros(n, 1);
F(1) = b(1) - b0;

k     = (2:n-1)';
rk    = r(k);
bPrev = b(k-1);
bCur  = b(k);
bNext = b(k+1);

F(k) = -lambda * ((bNext - 2 * bCur + bPrev) / h^2 ...
                + (2 ./ rk) .* (bNext - bPrev) / (2 * h)) ...
     + (bCur - 1) .* (nu * bCur + eta) ...
     + bCur ./ rk;

F(n) = b(n) - buffer.farField(r(n), lambda, nu, eta, 4);

lowerDiagonal = -lambda * (1 / h^2 - 1 ./ (rk * h));
mainDiagonal  =  2 * lambda / h^2 + 2 * nu * bCur + eta - nu + 1 ./ rk;
upperDiagonal = -lambda * (1 / h^2 + 1 ./ (rk * h));

J = spdiags([[lowerDiagonal; 0; 0], ...
             [1; mainDiagonal; 1], ...
             [0; 0; upperDiagonal]], -1:1, n, n);
end
