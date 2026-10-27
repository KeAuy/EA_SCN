function [residual, jacobian] = buffer_residual_jacobian(b, r, lambda, nu, eta, b0)
%BUFFER_RESIDUAL_JACOBIAN Residual and Jacobian of the finite-difference system.
%
%   The buffer equation is discretised on the uniform grid R with central
%   differences of second order, which gives the nonlinear system F(b) = 0
%   solved by SOLVE_BUFFER_BVP. This function returns F and its exact
%   Jacobian dF/db, which is tridiagonal.
%
%   Boundary conditions:
%       row 1:  b(r_min) = b0
%       row N:  b(r_max) = 1 - q/r + eta q^3/r^2
%                          + eta (1 - 2 q eta) q^4/r^3
%                          + eta (2 lambda + 5 eta q^2 (eta q - 1) q^4)/r^4
%   where q = 1/(eta + nu).
%
%   See also SOLVE_BUFFER_BVP.

q = 1 / (eta + nu);
n = numel(b);
h = r(2) - r(1);

residual = zeros(n, 1);

% Inner boundary: fixed value b0 at r_min.
residual(1) = b(1) - b0;

% Interior points: second-order central differences.
i     = (2:n-1).';
ri    = r(i);
bPrev = b(i-1);
bCur  = b(i);
bNext = b(i+1);

residual(i) = -lambda * ((bNext - 2 * bCur + bPrev) / h^2 ...
                         + (2 ./ ri) .* (bNext - bPrev) / (2 * h)) ...
              + (bCur - 1) .* (nu * bCur + eta) ...
              + bCur ./ ri;

% Outer boundary: truncated asymptotic expansion evaluated at r_max.
rMax  = r(end);
bInf  = 1 - q / rMax + eta * q^3 / rMax^2 ...
        + eta * (1 - 2 * q * eta) * q^4 / rMax^3 ...
        + eta * (2 * lambda + 5 * eta * q^2 * (eta * q - 1) * q^4) / rMax^4;
residual(n) = b(n) - bInf;

% Exact Jacobian of the interior stencil. The boundary rows are identity
% rows, so the matrix is tridiagonal.
lowerDiagonal = -lambda * (1 / h^2 - 1 ./ (ri * h));
mainDiagonal  =  2 * lambda / h^2 + 2 * nu * bCur + eta - nu + 1 ./ ri;
upperDiagonal = -lambda * (1 / h^2 + 1 ./ (ri * h));

% The three diagonals are padded to length n. SPDIAGS takes the subdiagonal
% from the first n-1 entries of its column and the superdiagonal from the
% last n-1 entries, hence the asymmetric padding.
jacobian = spdiags([[lowerDiagonal; 0; 0], ...
                    [1; mainDiagonal; 1], ...
                    [0; 0; upperDiagonal]], -1:1, n, n);
end
