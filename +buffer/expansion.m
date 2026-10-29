%EXPANSION Small-r (Frobenius) and large-r (far-field) series of the solution.
%
%   [BNEAR, DBNEAR, BFAR, DBFAR] = BUFFER.EXPANSION(R, B0, LAMBDA, NU, ETA,
%   NUMTERMS) returns NUMEL(R)-by-NUMEL(NU) arrays.  The near expansion is
%   b = b0 + b0/(2 lambda) r + beta r^2 with beta taken from the ODE.

function [bNear, dbNear, bFar, dbFar] = expansion(r, b0, lambda, nu, eta, numTerms)

arguments
    r double {mustBePositive}
    b0 double
    lambda (1,1) double {mustBePositive}
    nu double
    eta (1,1) double
    numTerms (1,1) double {mustBeMember(numTerms, [2 3 4])} = 4
end

r  = r(:);
b0 = b0(:)';
nu = nu(:)';

beta = ((b0 - 1) .* (nu .* b0 + eta) + b0 / (2 * lambda)) / (6 * lambda);
bNear  = b0 + (b0 / (2 * lambda)) .* r + beta .* r.^2;
dbNear = b0 / (2 * lambda) + 2 * beta .* r;

[bFar, dbFar] = buffer.farField(r, lambda, nu, eta, numTerms);
end
