function approx = approximationsGrid(r, lambdaValues, nuValues, eta, names, options)
%APPROXIMATIONSGRID Evaluate approximations on a lambda-nu parameter grid.
%
%   APPROX = BUFFER.APPROXIMATIONSGRID(R, LAMBDAVALUES, NUVALUES, ETA,
%   NAMES, OPTIONS) evaluates the approximations listed in NAMES for every
%   (lambda, nu) pair.  Every field of APPROX is NR-by-NLAMBDA-by-NNU.
%
%   The implementation loops over lambda and reuses BUFFER.APPROXIMATIONS,
%   so the closed-form expressions are maintained in one place.
%
%   OPTIONS are forwarded to BUFFER.APPROXIMATIONS: PadeRefine,
%   PadeRefineMethod, PadeRealCoefficients and RealBranchIndex.
%
%   See also BUFFER.APPROXIMATIONS.

arguments
    r double {mustBePositive}
    lambdaValues double
    nuValues double
    eta (1,1) double
    names = {}
    options.PadeRefine (1,1) logical = false
    options.PadeRefineMethod (1,:) char ...
        {mustBeMember(options.PadeRefineMethod, {'fsolve', 'newton'})} = 'fsolve'
    options.PadeRealCoefficients (1,1) logical = true
    options.RealBranchIndex (1,1) double {mustBeNonnegative} = 0
end

lambdaValues = lambdaValues(:);
nuValues     = nuValues(:);
nLambda      = numel(lambdaValues);
nNu          = numel(nuValues);
nR           = numel(r);

approx = struct();

for i = 1:nLambda
    current = buffer.approximations(r, lambdaValues(i), nuValues, eta, ...
        names, ...
        'PadeRefine', options.PadeRefine, ...
        'PadeRefineMethod', options.PadeRefineMethod, ...
        'PadeRealCoefficients', options.PadeRealCoefficients, ...
        'RealBranchIndex', options.RealBranchIndex);

    fields = fieldnames(current);
    for k = 1:numel(fields)
        field = fields{k};
        if i == 1
            approx.(field) = zeros(nR, nLambda, nNu);
        end
        approx.(field)(:, i, :) = reshape(current.(field), nR, 1, nNu);
    end
end
end
