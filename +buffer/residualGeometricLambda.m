%RESIDUALGEOMETRICLAMBDA Residual for every lambda column of B.
%
%   B is N-by-NLAMBDA and LAMBDA, B0 have one entry per column.  The layout
%   matches BUFFER.JACOBIANGEOMETRICLAMBDA and the batched Newton solver.

function F = residualGeometricLambda(b, r, hLeft, hRight, b0, lambda, nu, eta)

[n, m] = size(b);
lambda = lambda(:).';
b0     = b0(:).';
if numel(lambda) ~= m || numel(b0) ~= m
    error('buffer:sizeMismatch', ...
        'LAMBDA and B0 must have one entry per column of B.');
end

q = 1 / (eta + nu);

F = zeros(n, m);
F(1, :) = b(1, :) - b0;

k  = (2:n-1).';
hb = 0.5 * (hLeft(k) + hRight(k));

F(k, :) = lambda .* ((b(k+1, :) - b(k, :)) ./ hRight(k) ...
                   - (b(k, :) - b(k-1, :)) ./ hLeft(k)) ./ hb ...
        + 2 * lambda .* (b(k+1, :) - b(k-1, :)) ...
          ./ ((hLeft(k) + hRight(k)) .* r(k)) ...
        - (b(k, :) - 1) .* (nu * b(k, :) + eta) ...
        - b(k, :) ./ r(k);

F(n, :) = b(n, :) - (1 - q / r(n) + eta * q^3 / r(n)^2);
end
