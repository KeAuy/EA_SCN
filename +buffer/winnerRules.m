%WINNERRULES Method-region rules used by Figure 7.
%
%   WINNER = BUFFER.WINNERRULES(LAMBDAVALUES, NUVALUES, ETA, ALPHAGLOBAL,
%   OPTIONS) returns 1 = Pade2, 2 = RBA2, 3 = DblExp-Global.  The default is
%   Pade2; RBA2 is selected where lambda*nu < RBAProduct and lambda*eta <
%   RBALambda, and DblExp-Global where its alpha is real to within
%   REALALPHATOLERANCE.

function winner = winnerRules(lambdaValues, nuValues, eta, alphaGlobal, options)

arguments
    lambdaValues double
    nuValues double
    eta (1,1) double
    alphaGlobal double
    options.RBAProduct (1,1) double = 0.1
    options.RBALambda (1,1) double = 0.03
    options.RealAlphaTolerance (1,1) double {mustBePositive} = 1e-8
end

lambdaValues = lambdaValues(:);
nuValues     = nuValues(:);
nLambda      = numel(lambdaValues);
nNu          = numel(nuValues);

if ~isequal(size(alphaGlobal), [nLambda, nNu])
    error('buffer:sizeMismatch', ...
        'ALPHAGLOBAL must be NLAMBDA-by-NNU.');
end

winner = ones(nLambda, nNu);

isRBA = (lambdaValues * nuValues.' < options.RBAProduct) ...
        & (lambdaValues * eta < options.RBALambda);
winner(isRBA) = 2;

isRealGlobal = abs(imag(alphaGlobal)) < options.RealAlphaTolerance;
winner(~isRBA & isRealGlobal) = 3;
end
