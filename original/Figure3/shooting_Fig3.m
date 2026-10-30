clear; clc;
%%
% Fig2_A: lambda = 0.1, nu = 0.1
% Fig2_B: lambda = 0.1, nu = 10
% Fig2_C: lambda = 1, nu = 10
%%
lambda = 1;
nu = 10;
eta = 1;
N = 100;
b0_test = 0.5;

n = 1:N;
r = 10.^(-3 + 5*n/N);

options = bvpset('NMax', 50000, 'RelTol', 1e-10, 'AbsTol', 1e-10);
yinit   = bvpinit(r, @(rr) guess(b0_test, lambda, nu, eta, rr));

sol = bvp5c(@(rr, y) odefun(rr, y, nu, lambda, eta), ...
            @(ya, yb) bcfun(ya, yb, lambda), ...
            yinit, options);

Y = deval(sol, r);
fprintf('nu = %.3e completed\n', nu);
fprintf('initial value b0: %.3e\n', Y(1,1));
