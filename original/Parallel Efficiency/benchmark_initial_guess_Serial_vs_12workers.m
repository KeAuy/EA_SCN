clear; clc; close all;
%% Configuration

nLambda = 100;
nNu     = 100;
nGrid   = 100;
eta     = 1;
tol     = 1e-6;
nMax    = 50000;

lambdaValues = logspace(-3, 2, nLambda)';
nuValues     = logspace(-3, 2, nNu)';
workersList  = [1 2 4 8 12];
nRepeats     = 3;

%%
outDir = fileparts(mfilename('fullpath'));
cacheFile  = fullfile(outDir, 'benchmark_b0_100x100.mat');

if isfile(cacheFile)
    S = load(cacheFile);
    cacheMatches = isfield(S, 'results') && ~isempty(S.results) ...
        && isfield(S, 'nLambda') && S.nLambda == nLambda ...
        && isfield(S, 'nNu')     && S.nNu == nNu ...
        && isfield(S, 'nGrid')   && S.nGrid == nGrid ...
        && isfield(S, 'eta')     && S.eta == eta ...
        && isfield(S, 'tol')     && S.tol == tol ...
        && numel(S.results) == numel(workersList) ...
        && isequal([S.results.workers], workersList);

    if cacheMatches
        fprintf('Loaded cached results from: %s\n', cacheFile);
        reportGenerationResults(S.results);
        return
    end
end
