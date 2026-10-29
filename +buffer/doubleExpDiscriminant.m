%DOUBLEEXPDISCRIMINANT Cubic discriminant of the double-exponential alpha.
%
%   DELTA = BUFFER.DOUBLEEXPDISCRIMINANT(LAMBDAVALUES, NUVALUES, ETA, TYPE)
%   returns the normalized discriminant; its zero contour separates the real
%   and complex alpha branches and is drawn on the Figure 6/7 method maps.

function delta = doubleExpDiscriminant(lambdaValues, nuValues, eta, type)

arguments
    lambdaValues double
    nuValues double
    eta (1,1) double
    type (1,:) char {mustBeMember(type, {'var', 'global'})}
end

lambdaValues = lambdaValues(:);
nuValues     = nuValues(:);
nLambda      = numel(lambdaValues);
nNu          = numel(nuValues);
delta        = zeros(nLambda, nNu);

for i = 1:nLambda
    [pCoef, qCoef, rCoef, q] = buffer.doubleExpCoefficients( ...
        lambdaValues(i), nuValues, eta, type);

    a = eta * q.^2 .* pCoef;
    b = -qCoef;
    c = -rCoef;
    d = 1 ./ q;

    delta(i, :) = (18 * a .* b .* c .* d ...
                   - 4 * b.^3 .* d ...
                   + b.^2 .* c.^2 ...
                   - 4 * a .* c.^3 ...
                   - 27 * a.^2 .* d.^2).';
end

scale = max(abs(delta(:)));
if scale > 0
    delta = delta / scale;
end
end
