clear; clc;

lambda = 20.0;
eta = 1;
N = 100;

b0_test = 0.5;
nu_vec = logspace(-3, 2, N);
n = 1:N;
r = 10.^(-3 + 5*n/N);

options = bvpset('NMax', 50000, 'RelTol', 1e-8, 'AbsTol', 1e-8);
b0 = zeros(1, N);

nu = nu_vec(1);
sol = bvp5c(@(rr, y) odefun(rr, y, nu, lambda, eta), ...
            @(ya, yb) bcfun(ya, yb, lambda), ...
            bvpinit(r, @(rr) guess(b0_test, lambda, nu, eta, rr)), ...
            options);
Y = deval(sol, r);
b0(1) = Y(1,1);
fprintf('nu = %.3e completed\n', nu);

for k = 2:N
    nu = nu_vec(k);
    sol = bvp5c(@(rr, y) odefun(rr, y, nu, lambda, eta), ...
                @(ya, yb) bcfun(ya, yb, lambda), ...
                sol, options);
    Y = deval(sol, r);
    b0(k) = Y(1,1);
    fprintf('nu = %.3e completed\n', nu);
end

save('b0_Fig2_C.mat', 'b0');
