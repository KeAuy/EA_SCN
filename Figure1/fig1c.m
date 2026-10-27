%% Figure 1C: numerical solution and closed-form approximations of b(r)
%
% Same comparison as Figure 1A for lambda = 1, nu = 10, eta = 1, including
% the solution of the linearised equation.
%
% See README.md for the model and for the list of approximations.

clear; close all; clc;

% Model parameters
lambda = 1;
nu     = 10;
eta    = 1;
b0     = 0.735;   % b(r_min), fixed by the regularity condition at r = 0

% Numerical grid
rMin  = 1e-3;
rMax  = 1e2;
nGrid = 50000;
r     = linspace(rMin, rMax, nGrid)';

% Numerical solution
b = solve_buffer_bvp(b0, lambda, nu, eta, r, ...
    'SwitchRadius', 1.0, 'Tolerance', 1e-10, 'MaxIterations', 25);

% Closed-form approximations
approx = analytic_approximations(r, lambda, nu, eta);

% Plot style
gray  = [0.75 0.75 0.75];
green = [46 139 87] / 256;

figure;
hold on;
box off;

plot(r, b,                              'Color', gray,  'LineWidth', 2.5);
plot(r, approx.linearized,              '-.', 'Color', 'k',   'LineWidth', 1.25);
plot(r, approx.pade,                    '-',  'Color', green, 'LineWidth', 2.5);
plot(r, real(approx.pade2),             '--', 'Color', green, 'LineWidth', 2.5);
plot(r, approx.exp_series,              '-',  'Color', 'k',   'LineWidth', 2.5);
plot(r, approx.double_exp_series,       '--', 'Color', 'm',   'LineWidth', 2.5);

set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]);
xticks([1e-3 1e-2 1e-1 1 10 100]);
xlabel('r (distance)');
ylabel('b(r)');
title('Numerical solution and analytical approximations');
subtitle('\lambda = 1, \nu = 10, \eta = 1');
ylim([0.73, 1]);
yticks([0.73 0.75 0.8 0.85 0.9 0.95 1]);
yticklabels({' ', '0.75', '0.8', '0.85', '0.9', '0.95', '1'});
legend('Numerical', 'LIN', 'Padé', 'Padé2', 'Exp-Ser', 'DblExp-Ser', ...
       'Location', 'southeast');
