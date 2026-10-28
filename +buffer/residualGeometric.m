function F = residualGeometric(b, r, hLeft, hRight, b0, lambda, nu, eta)
%RESIDUALGEOMETRIC Residual of the finite-difference system on a geometric grid.
%
%   F = BUFFER.RESIDUALGEOMETRIC(B, R, HLEFT, HRIGHT, B0, LAMBDA, NU, ETA)
%   evaluates the second-order central-difference residual of
%
%       lambda (b'' + 2 b'/r) = (b - 1)(nu b + eta) + b/r
%
%   on the non-uniform grid R.  HLEFT(k) = R(k) - R(k-1) and
%   HRIGHT(k) = R(k+1) - R(k).  The inner value is B(1) = B0 and the outer
%   condition uses the two-term far-field expansion.
%
%   See also BUFFER.JACOBIANGEOMETRIC, BUFFER.SOLVEGEOMETRIC.

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
