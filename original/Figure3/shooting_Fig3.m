clear; clc;
%%
% Fig3_A: lambda = 0。1, nu = 0.1,eta = 1
% Fig3_B: lambda = 0.1, nu = 10, eta = 1
% Fig3_C: lambda = 1, nu = 10, eta = 1
%%
% parameters
lambda = 1; 
nu = 10; 
eta = 1;
N = 100;
b0_test = 0.5;

n = 1:100;
r = 10.^(-3+5*n/N);

% storage
sol_all = cell(N,1);
y_num = zeros(N,length(r));

options = bvpset('NMax',50000,'RelTol',1e-10,'AbsTol',1e-10);

yinit_struct = bvpinit(r, @(r)guess(b0_test, lambda, nu, eta, r));


sol = bvp5c(@(r,y) odefun(r, y, nu, lambda, eta), ...
            @(ya,yb) bcfun(ya, yb, lambda), ...
            yinit_struct);

sol_all{1} = sol;
sol_prev = sol;
Y = deval(sol,r);
b0(1) = Y(1,1);
y_num(:,1) = Y(1,:);
fprintf('nu = %.3e completed\n', nu);
fprintf('initial value b0: %.3e\n', b0(1));