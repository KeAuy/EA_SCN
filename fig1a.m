%% Figure 1A: numerical solution and closed-form approximations of b(r)
%
% Solves the buffer equation with the finite-difference/Newton solver and
% compares the result with the approximations used in the paper
% (lambda = 0.1, nu = 0.1, eta = 1).
%
% See README.md for the model and for the list of approximations.

clear; close all; clc;

% Model parameters
lambda = 0.1;
nu     = 0.1;
eta    = 1;
b0     = 0.131;   % b(r_min), fixed by the regularity condition at r = 0

% Numerical grid
rMin  = 1e-3;
rMax  = 1e2;
nGrid = 50000;
r     = linspace(rMin, rMax, nGrid)';

% Numerical solution
b = solve_buffer_bvp(b0, lambda, nu, eta, r, ...
    'SwitchRadius', 3.0, 'Tolerance', 1e-8, 'MaxIterations', 125);

% Closed-form approximations
approx = analytic_approximations(r, lambda, nu, eta);

% Plot style
gray  = [0.75 0.75 0.75];
green = [46 139 87] / 256;

figure;
hold on;
box off;

plot(r, b,                              'Color', gray,  'LineWidth', 2);
plot(r, approx.rba2,                    '--', 'Color', 'r',   'LineWidth', 2);
plot(r, approx.pade,                    '-',  'Color', green, 'LineWidth', 2);
plot(r, real(approx.pade2),             '--', 'Color', green, 'LineWidth', 2);
plot(r, approx.exp_series,              '-',  'Color', 'k',   'LineWidth', 2);
plot(r, real(approx.exp_pade),          '--', 'Color', 'k',   'LineWidth', 2);
plot(r, real(approx.double_exp_series), '--', 'Color', 'm',   'LineWidth', 2);

xlabel('r (distance)');
ylabel('b(r)');
title('Numerical solution and analytical approximations');
subtitle('\lambda = 0.1, \nu = 0.1, \eta = 1');
legend('Numerical', 'RBA2', 'Padé', 'Padé2', 'Exp-Ser', 'EXP-Padé', 'DblExp-Ser', ...
       'Location', 'northwest');

set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]);
xticks([1e-3 1e-2 1e-1 1 10 100]);
ylim([0.089, 1]);
yticks([0 0.2 0.4 0.6 0.8 1]);
yticklabels({' ', '0.2', '0.4', '0.6', '0.8', '1'});
