clear; clc;

eta     = 1;
b0_test = 0.5;
N     = 100;
N1      = 100;

r        = 10.^(-3 + 5*(1:N)/N);
lambda_v = logspace(-3, 2, N1);
nu_v     = logspace(-3, 2, N1);

options = bvpset('NMax', 50000, 'RelTol', 1e-5, 'AbsTol', 1e-5);

b0      = zeros(N1, N1);
b0p     = zeros(N1, N1);
y_num   = zeros(N1, N1, N);
sol_all = cell(N, N);

for i = 1:N1
    lam = lambda_v(i);
    for k = 1:N1                       % ← 这里循环 nu，不是 r
        nu = nu_v(k);

        if k > 1
            solinit = sol_all{i, k-1};        % 同一 λ 行，沿 ν 热启动
        elseif i > 1
            solinit = sol_all{i-1, 1};        % 上一行的起点
        else
            solinit = bvpinit(r, @(rr) guess(b0_test, lam, nu, eta, rr));
        end

        solij = bvp5c(@(rr,y) odefun(rr, y, nu, lam, eta), ...
                      @(ya,yb) bcfun(ya, yb, lam), ...
                      solinit, options);

        sol_all{i,k} = solij;
        Yij          = deval(solij, r);
        b0(i,k)      = Yij(1,1);
        b0p(i,k)     = Yij(2,1);
        y_num(i,k,:) = Yij(1,:);

        fprintf('lambda = %.3e, nu = %.3e 完成\n', lam, nu);
    end
end
save('b0_Fig6.mat','b0');