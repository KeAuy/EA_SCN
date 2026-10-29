%INITIALPROFILELAMBDA Initial profile for several lambda values at once.
%
%   B is NUMEL(R)-by-NUMLAMBDA; each column is built by BUFFER.INITIALPROFILE
%   for its own lambda with the same nu and eta.

function b = initialProfileLambda(r, b0, lambdaValues, nu, eta)
arguments
    r double {mustBePositive}
    b0 double
    lambdaValues double
    nu (1,1) double
    eta (1,1) double
end

r            = r(:);
b0           = b0(:).';
lambdaValues = lambdaValues(:);
nLambda      = numel(lambdaValues);

if numel(b0) ~= nLambda
    error('buffer:sizeMismatch', ...
        'B0 must have one entry per LAMBDA.');
end

b = zeros(numel(r), nLambda);

for i = 1:nLambda
    b(:, i) = buffer.initialProfile(r, b0(i), lambdaValues(i), nu, eta, []);
end
end
