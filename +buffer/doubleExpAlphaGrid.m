function alpha = doubleExpAlphaGrid(lambdaValues, nuValues, eta, type)
%DOUBLEEXPALPHAGRID Double-exponential alpha on a lambda-nu grid.
%
%   ALPHA = BUFFER.DOUBLEEXPALPHAGRID(LAMBDAVALUES, NUVALUES, ETA, TYPE)
%   returns an NLAMBDA-by-NNU matrix with the first cubic root used by the
%   Exp-Var or Exp-Global approximation.  TYPE is 'var' or 'global'.
%
%   See also BUFFER.DOUBLEEXPCOEFFICIENTS, BUFFER.WINNERRULES.

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
