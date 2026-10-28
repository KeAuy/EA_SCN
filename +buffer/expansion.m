function [bNear, dbNear, bFar, dbFar] = expansion(r, b0, lambda, nu, eta, numTerms)
%EXPANSION Small-r and far-field series expansions of the buffer solution.
%
%   [BNEAR, DBNEAR, BFAR, DBFAR] = BUFFER.EXPANSION(R, B0, LAMBDA, NU, ETA,
%   NUMTERMS) evaluates the Frobenius expansion near r = 0 and the truncated
%   far-field expansion of
%
%       b'' + (2/r) b' - (1/lambda) (b - 1) (nu b + eta) - b/(lambda r) = 0
%
%   together with their first derivatives.  R is a vector; B0 and NU may be
%   scalars or vectors.  The outputs are NUMEL(R)-by-NUMEL(NU).
%
%   NUMTERMS selects the number of far-field terms (2, 3 or 4).
%
%   See also BUFFER.FARFIELD, BUFFER.GUESS, BUFFER.INITIALPROFILE.

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

% Small-r Frobenius expansion: b = b0 + b0/(2 lambda) r + beta r^2 + O(r^3).
beta = ((b0 - 1) .* (nu .* b0 + eta) + b0 / (2 * lambda)) / (6 * lambda);
bNear  = b0 + (b0 / (2 * lambda)) .* r + beta .* r.^2;
dbNear = b0 / (2 * lambda) + 2 * beta .* r;

[bFar, dbFar] = buffer.farField(r, lambda, nu, eta, numTerms);
end
