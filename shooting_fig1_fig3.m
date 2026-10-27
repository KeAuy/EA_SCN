%% Shooting solution of the buffer equation with BVP5C
%
% Solves the boundary value problem
%
%   b'' + (2/r) b' - (1/lambda) (b - 1) (nu b + eta) - b / (lambda r) = 0,
%   b'(r_min) = b(r_min) / (2 lambda),   b(r_max) = 1
%
% on a logarithmic grid, using the series expansions of BUFFER_GUESS as the
% initial guess. This is the reference solution used by the shooting
% continuation scripts.
%
% See README.md for the model and for the shooting interface
% (BUFFER_ODE, BUFFER_BC, BUFFER_GUESS).

clear; clc;

% Model parameters
lambda  = 1;
nu      = 10;
eta     = 1;
b0_test = 0.5;   % parameter of the initial guess

% Logarithmic grid from 1e-3 to 1e2
nGrid = 100;
r = 10.^(-3 + 5 * (1:nGrid) / nGrid);

% Solve the boundary value problem
options      = bvpset('NMax', 50000, 'RelTol', 1e-10, 'AbsTol', 1e-10);
initialGuess = bvpinit(r, @(rr) buffer_guess(rr, b0_test, lambda, nu, eta));
solution     = bvp5c(@(rr, y) buffer_ode(rr, y, lambda, nu, eta), ...
                     @(ya, yb) buffer_bc(ya, yb, lambda), ...
                     initialGuess, options);

% Value at the inner boundary
y = deval(solution, r);
fprintf('nu = %.3e completed\n', nu);
fprintf('b(r_min) = %.3e\n', y(1, 1));
