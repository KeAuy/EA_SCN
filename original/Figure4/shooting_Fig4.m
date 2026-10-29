clear; clc;
%%
% Fig4_A: lambda = 0.02,eta = 1
% Fig4_B: lambda = 2, eta = 1
% Fig4_C: lambda = 20, eta = 1,N = 1000
%%
lambda = 20;
eta = 1;
N = 1000;

b0_test = 0.5;
nu_vec = logspace(-3, 2, N);
nu = nu_vec(1);
n = 1:N;
r = 10.^(-3+5*n/N);

% storage
sol_all = cell(N,1);
b0 = zeros(1,N);
y_num = zeros(N,length(r));

options = bvpset('NMax',50000,'RelTol',1e-8,'AbsTol',1e-8);

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

for k = 2:N
    nu = nu_vec(k);

    solinit = sol_prev;
    sol = bvp5c(@(r,y) odefun(r,y,nu,lambda,eta), ...
                @(ya,yb) bcfun(ya, yb, lambda), ...
                solinit, ...
                options);

    sol_all{k} = sol;
    y = deval(sol, r);
    b0(k) = y(1,1);
    y_num(:,k) = y(1,:);
    sol_prev = sol;
    fprintf('nu = %.3e completed\n', nu);
end
save('b0_Fig4_C.mat',"b0");