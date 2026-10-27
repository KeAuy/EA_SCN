function [bInner, dbInner, bOuter, dbOuter] = buffer_expansion(r, b0, lambda, nu, eta)
%BUFFER_EXPANSION Series expansions of the buffer solution b(r).
%
%   [BINNER, DBINNER, BOUTER, DBOUTER] = BUFFER_EXPANSION(R, B0, LAMBDA, NU, ETA)
%   evaluates the small-r (Frobenius) and large-r (asymptotic) expansions of
%   the solution of the buffer equation
%
%       b'' + (2/r) b' - (1/lambda) (b - 1) (nu b + eta) - b / (lambda r) = 0
%
%   together with their first derivatives db/dr.
%
%   The two branches are combined by BUFFER_GUESS (initial guess for BVP5C)
%   and by SOLVE_BUFFER_BVP (initial profile for the Newton iteration).
%
%   See also BUFFER_GUESS, SOLVE_BUFFER_BVP.

q    = 1 / (eta + nu);
beta = ((b0 - 1) * (nu * b0 + eta) + b0 / (2 * lambda)) / (6 * lambda);

% Small-r Frobenius expansion (regular singular point at r = 0):
%   b = b0 + b0/(2 lambda) r + beta r^2 + O(r^3)
bInner  = b0 + (b0 / (2 * lambda)) * r + beta * r.^2;
dbInner = b0 / (2 * lambda) + 2 * beta * r;

% Large-r asymptotic expansion in powers of x = 1/r:
%   b = 1 - q x + eta q^3 x^2 + eta (1 - 2 q eta) q^4 x^3
%       + eta (2 lambda + 5 eta q^2 (eta q - 1)) q^4 x^4 + O(x^5)
x      = 1 ./ r;
bOuter = 1 - q * x + eta * q^3 * x.^2 ...
         + eta * (1 - 2 * q * eta) * q^4 * x.^3 ...
         + eta * (2 * lambda + 5 * eta * q^2 * (eta * q - 1)) * q^4 * x.^4;

% db/dr = (db/dx) (dx/dr) with dx/dr = -1/r^2.
dbdx    = -q + 2 * eta * q^3 * x + 3 * eta * (1 - 2 * q * eta) * q^4 * x.^2 ...
          + 4 * eta * (2 * lambda + 5 * eta * q^2 * (eta * q - 1)) * q^4 * x.^3;
dbOuter = -dbdx ./ r.^2;
end
