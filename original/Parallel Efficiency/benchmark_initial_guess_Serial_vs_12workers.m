clear; clc; close all;

%%  Configuration
nLambda  = 100;
nNu      = 100;
nGrid    = 100;
eta      = 1;
tol      = 1e-6;
nMax     = 50000;
nWorkers = 12;
nRepeats = 3;

lambdaValues = logspace(-3, 2, nLambda)';
nuValues     = logspace(-3, 2, nNu)';


%%
outDir = fileparts(mfilename('fullpath'));
cacheFile  = fullfile(outDir, 'benchmark_b0_serial_vs_12.mat');

if isfile(cacheFile)
    S = load(cacheFile);
    cacheMatches = isfield(S, 'result') ...
        && isfield(S, 'nLambda') && S.nLambda == nLambda ...
        && isfield(S, 'nNu')     && S.nNu == nNu ...
        && isfield(S, 'nGrid')   && S.nGrid == nGrid ...
        && isfield(S, 'eta')     && S.eta == eta ...
        && isfield(S.result, 'workers') && S.result.workers == nWorkers;

    if cacheMatches
        fprintf('Loaded cached results from: %s\n', cacheFile);
        report_Serial_vs_12workers(S.result);
        return
    end
end
