%% Figure 6 B1/B2: best-method map and c-error map, eta = 1.
%
% Same parameter grid as Figure 6 A1/A2 but with a 1400-point radial grid and
% the charge-density error c = nu (b - 1) + 1/r.  Exp-Pade is masked out for
% nu >= eta.
%
% Output: figures/figure6B1B2_winnerMap.png
%         figures/figure6B1B2_cErrorMap.png

clear; close all; clc;

figureDir = fileparts(mfilename('fullpath'));
addpath(fullfile(figureDir, '..'));
outDir = buffer.outputDir();

eta           = 1;
nGrid         = 1400;
nLambda       = 100;
nNu           = 100;
tolerance     = 1e-3;
stepTolerance = 1e-12;
maxIterations = 500;

lambdaValues = logspace(-3, 2, nLambda)';
nuValues     = logspace(-3, 2, nNu)';
grid         = buffer.geometricGrid(1e-3, 1e2, nGrid);
b0           = buffer.generateB0Grid( ...
    lambdaValues, nuValues, eta, 'Verbose', false);

bNum = buffer.solveGeometricLambdaGrid(b0, lambdaValues, nuValues, ...
    eta, grid, 'Tolerance', tolerance, 'StepTolerance', stepTolerance, ...
    'MaxIterations', maxIterations, 'Verbose', false);

methodNames = {'Pade2', 'RBA2', 'Exp-Pade', 'Dbl-ExpGlobal', 'DblExp-Var'};
approx = buffer.approximationsGrid(grid.r, lambdaValues, nuValues, eta, ...
    {'Pade2', 'RBA2', 'ExpPade', 'DblExpGlobal', 'DblExpVar'}, ...
    'PadeRefine', true, 'PadeRefineMethod', 'newton', ...
    'PadeRealCoefficients', true);

cNum = buffer.chargeDensityGrid(bNum, grid.r, nuValues);
cApprox = { ...
    buffer.chargeDensityGrid(approx.pade2, grid.r, nuValues), ...
    buffer.chargeDensityGrid(approx.rba2, grid.r, nuValues), ...
    buffer.chargeDensityGrid(approx.expPade, grid.r, nuValues), ...
    buffer.chargeDensityGrid(approx.doubleExpGlobal, grid.r, nuValues), ...
    buffer.chargeDensityGrid(approx.doubleExpVar, grid.r, nuValues)};

badNu = find(nuValues >= eta);
if ~isempty(badNu)
    cApprox{3}(:, :, badNu) = NaN;
end

errors = cell(1, numel(cApprox));
for k = 1:numel(cApprox)
    value    = real(cApprox{k});
    valid    = (value > 0) & (cNum > 0);
    errors{k} = buffer.meanErrorGrid(cNum, value, 'log', ...
        'ValidMask', valid, 'InvalidValue', Inf);
end
errors = cat(3, errors{:});

[winner, bestError] = buffer.winnerMap(errors);

colors = [0.13 0.69 0.30; ...
          0.55 0.10 0.10; ...
          0.50 0.50 0.50; ...
          1.00 0.61 0.82; ...
          0.90 0.00 0.90];

deltaVar    = buffer.doubleExpDiscriminant(lambdaValues, nuValues, eta, 'var');
deltaGlobal = buffer.doubleExpDiscriminant(lambdaValues, nuValues, eta, 'global');
contours = struct( ...
    'Z',         {deltaVar, deltaGlobal}, ...
    'Color',     {[0.85 0.85 0.85], [0.85 0.65 0.00]}, ...
    'LineWidth', {1.0, 1.0});

buffer.plotWinnerMap(nuValues, lambdaValues, winner, methodNames, colors, ...
    'Title', sprintf('Method regions (\\eta = %g)', eta), ...
    'Contours', contours, ...
    'OutputFile', fullfile(outDir, 'figure6B1B2_winnerMap.png'));

errorMap = log10(bestError);
errorMap = min(max(errorMap, -9), -2);

buffer.plotErrorMap(nuValues, lambdaValues, errorMap, ...
    'Title', 'log_{10}||c_{num} - c_{approx}||, \eta = 1', ...
    'CLim', [-9 -2], 'ColorTicks', -9:1:-2, ...
    'OutputFile', fullfile(outDir, 'figure6B1B2_cErrorMap.png'));
