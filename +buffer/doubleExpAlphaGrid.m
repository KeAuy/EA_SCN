%DOUBLEEXPALPHAGRID Double-exponential alpha on a lambda-nu grid.
%
%   ALPHA = BUFFER.DOUBLEEXPALPHAGRID(LAMBDAVALUES, NUVALUES, ETA, TYPE)
%   returns an NLAMBDA-by-NNU matrix of the first cubic root for TYPE 'var'
%   or 'global'.

function alpha = doubleExpAlphaGrid(lambdaValues, nuValues, eta, type)

arguments
    lambdaValues double
    nuValues double
    eta (1,1) double
    type (1,:) char {mustBeMember(type, {'var', 'global'})}
end

lambdaValues = lambdaValues(:);
nuValues     = nuValues(:);
nLambda      = numel(lambdaValues);
alpha        = zeros(nLambda, numel(nuValues));

for i = 1:nLambda
    [pCoef, qCoef, rCoef, q] = buffer.doubleExpCoefficients( ...
        lambdaValues(i), nuValues, eta, type);
    alpha(i, :) = buffer.doubleExpAlpha(q, eta, pCoef, qCoef, rCoef).';
end
end
