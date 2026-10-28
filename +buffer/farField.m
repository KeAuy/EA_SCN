function [b, db] = farField(r, lambda, nu, eta, numTerms)
%FARFIELD Truncated large-r expansion of the buffer solution and its derivative.
%
%   [B, DB] = BUFFER.FARFIELD(R, LAMBDA, NU, ETA, NUMTERMS) evaluates
%
%       b = 1 - q x + eta q^3 x^2
%           + eta (1 - 2 eta q) q^4 x^3
%           + eta (2 lambda + 5 eta q^2 (eta q - 1)) q^4 x^4
%
%   with x = 1/r and q = 1/(nu + eta), using NUMTERMS = 2, 3 or 4 terms.
%
%   See also BUFFER.EXPANSION.

arguments
    r double {mustBePositive}
    lambda (1,1) double {mustBePositive}
    nu double
    eta (1,1) double
    numTerms (1,1) double {mustBeMember(numTerms, [2 3 4])}
end

r  = r(:);
nu = nu(:)';
q  = 1 ./ (eta + nu);
x  = 1 ./ r;

b   = 1 - q .* x + eta * q.^3 .* x.^2;
dbdx = -q + 2 * eta * q.^3 .* x;

if numTerms >= 3
    b    = b + eta * (1 - 2 * eta * q) .* q.^4 .* x.^3;
    dbdx = dbdx + 3 * eta * (1 - 2 * eta * q) .* q.^4 .* x.^2;
end

if numTerms >= 4
    c4   = eta * (2 * lambda + 5 * eta * q.^2 .* (eta * q - 1)) .* q.^4;
    b    = b + c4 .* x.^4;
    dbdx = dbdx + 4 * c4 .* x.^3;
end

% Chain rule: db/dr = (db/dx) (dx/dr) with dx/dr = -1/r^2.
db = -dbdx ./ r.^2;
end
