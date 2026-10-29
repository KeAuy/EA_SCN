%APPROXIMATIONSGRID Evaluate the approximations on a lambda-nu grid.
%
%   APPROX = BUFFER.APPROXIMATIONSGRID(R, LAMBDAVALUES, NUVALUES, ETA, NAMES,
%   OPTIONS) loops over lambda and stores every field as NR-by-NLAMBDA-by-NNU.
%   The closed-form expressions themselves live in BUFFER.APPROXIMATIONS.

function approx = approximationsGrid(r, lambdaValues, nuValues, eta, names, options)

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
