%GENERATEB0GRID b0(lambda, nu) on a parameter grid.
%
%   B0 = BUFFER.GENERATEB0GRID(LAMBDAVALUES, NUVALUES, ETA, OPTIONS) solves
%   the buffer BVP for every lambda, continuing along nu inside each row.  The
%   result is NLAMBDA-by-NNU and feeds the lambda-batch solver of Figures 6/7.
%
%   With the Parallel Computing Toolbox the rows run in a PARFOR loop, one
%   worker per lambda, and a progress bar is drawn from a DataQueue.  Without
%   the toolbox, or for grids with fewer than eight lambdas, the loop stays
%   serial.
%
%   OPTIONS: GridSize, Tolerance, NMax, InitialB0, UseParallel, Progress and
%   Verbose.  Nothing is written to disk.

function b0 = generateB0Grid(lambdaValues, nuValues, eta, options)

arguments
    lambdaValues double
    nuValues double
    eta (1,1) double
    options.GridSize = []
    options.Tolerance (1,1) double {mustBePositive} = 1e-6
    options.NMax (1,1) double {mustBeInteger, mustBePositive} = 50000
    options.InitialB0 (1,1) double = 0.5
    options.UseParallel (1,1) logical = true
    options.Progress (1,1) logical = true
    options.Verbose (1,1) logical = true
end

lambdaValues = lambdaValues(:);
nuValues     = nuValues(:);
nLambda      = numel(lambdaValues);
nNu          = numel(nuValues);

if isempty(options.GridSize)
    options.GridSize = nNu;
end
if options.GridSize < 3
    error('buffer:invalidGrid', 'GridSize must be at least three.');
end

grid      = buffer.geometricGrid(1e-3, 1e2, options.GridSize);
r         = grid.r;
tolerance = options.Tolerance;
nMax      = options.NMax;
initialB0 = options.InitialB0;
b0        = zeros(nLambda, nNu);

useParallel = options.UseParallel && nLambda >= 8 && hasParallelToolbox();
if useParallel
    try
        if isempty(gcp('nocreate'))
            parpool('local');
        end
    catch
        useParallel = false;
    end
end

if useParallel
    queue = parallel.pool.DataQueue;
    if options.Progress
        afterEach(queue, buffer.makeProgress(nLambda, lambdaValues, tic));
    end

    parfor i = 1:nLambda
        b0(i, :) = buffer.generateB0Row(lambdaValues(i), nuValues, eta, ...
            r, tolerance, nMax, initialB0);
        if options.Progress
            send(queue, i);
        end
    end
else
    for i = 1:nLambda
        b0(i, :) = buffer.generateB0Row(lambdaValues(i), nuValues, eta, ...
            r, tolerance, nMax, initialB0);
        if options.Verbose
            fprintf('lambda = %.3e completed\n', lambdaValues(i));
        end
    end
end
end

function tf = hasParallelToolbox()
tf = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
end
