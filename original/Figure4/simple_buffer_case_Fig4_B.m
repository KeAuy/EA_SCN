clear; close all; clc;

lambda   = 2;
eta      = 1.0;
N        = 100;

tol      = 1e-7;
max_iter = 1000;
nu_vec = logspace(-3, 2, N)';

k_idx  = (1:N)';
r_vec  = 10.^(-3 + 5.*k_idx./N);
q_vec  = 1 ./ (eta + nu_vec);
x_vec  = 1 ./ r_vec;

alpha_r = 10^(5/N);
h_L = r_vec .* (1 - 1/alpha_r);
h_R = r_vec .* (alpha_r - 1);

temp  = load("b0_Fig4_B.mat");
b_r0  = temp.b0;

b_frob = zeros(N, N);

for j = 1:N
    b0j   = b_r0(j);
    nu_j  = nu_vec(j);
    c1_j  = b0j / (2*lambda);
    c2_j  = ((b0j - 1)*(nu_j*b0j + eta) + c1_j) / (6*lambda);
    b_frob(:, j) = b0j + c1_j .* r_vec + c2_j .* (r_vec.^2);
end

b_asymp = zeros(N, N);

for j = 1:N
    q_j = q_vec(j);
    b_asymp(:, j) = 1 ...
        - q_j   .* x_vec ...
        + eta*(q_j^3)         .* x_vec.^2 ...
        + eta*(1 - 2*eta*q_j)*(q_j^4) .* x_vec.^3 ...
        + eta*(2*lambda + 5*eta*(q_j^2)*(eta*q_j - 1))*(q_j^4) .* x_vec.^4;
end

rel_disc   = abs(b_frob - b_asymp) ./ (abs(b_frob) + abs(b_asymp));
[~, idx_switch] = min(rel_disc, [], 1);
r_switch   = r_vec(idx_switch);

b_init = zeros(N, N);

for j = 1:N
    mask = r_vec < r_switch(j);

    b_init(:, j) = double(mask).*b_frob(:, j) ...
                 + double(~mask).*b_asymp(:, j);
end

b_num       = b_init;
res_vec     = zeros(N, 1);
failed_cols = zeros(2, N);

for j = 1:N
    b0j        = b_r0(j);
    b_col_init = b_init(:, j);
    converged  = false;

    for iter = 1:max_iter
        F_col = compute_F(b0j, b_num(:,j), r_vec, h_L, h_R, N, lambda, nu_vec(j), eta);
        res_vec(j) = norm(F_col, inf);
        if res_vec(j) < tol
            fprintf('  Newton 收敛，j=%d\n', j);
            converged = true;
            break;
        end
        J_col = compute_J(b_num(:,j), r_vec, h_L, h_R, N, lambda, nu_vec(j), eta);
        delta = -J_col \ F_col;
        b_num(:, j) = b_num(:, j) + delta;
    end

    if ~converged
        fprintf('  Newton False，j=%d，引入松弛项。\n', j);
        b_num(:, j) = b_col_init;

        for iter = 1:max_iter
            F_col   = compute_F_Fig2(b0j, b_num(:,j), r_vec, h_L, h_R, N, lambda, nu_vec(j), eta);
            res_old = norm(F_col, inf);
            if res_old < tol
                fprintf(' 松弛法迭代结果收敛，j=%d\n', j);
                converged = true;
                break;
            end
            J_col = compute_J(b_num(:,j), r_vec, h_L, h_R, N, lambda, nu_vec(j), eta);
            delta = -J_col \ F_col;

            omega_grid = 0.1:0.05:1.0;
            best_res   = inf;
            best_b     = b_num(:, j);
            for omega = omega_grid
                b_trial   = b_num(:, j) + omega * delta;
                F_trial   = compute_F_Fig2(b0j, b_trial, r_vec, h_L, h_R, N, lambda, nu_vec(j), eta);
                res_trial = norm(F_trial, inf);
                if res_trial < best_res
                    best_res = res_trial;
                    best_b   = b_trial;
                end
            end

            if best_res < res_old
                b_num(:, j) = best_b;
                res_vec(j)  = best_res;
            else
                fprintf('松弛法失败，j=%d, iter=%d\n', j, iter);
                failed_cols(1, j) = j;
                failed_cols(2, j) = nu_vec(j);
                break;
            end
        end
    end

    if ~converged
        warning('第 %d 列 未收敛', j);
    end
end
save("b_num_Fig4_B.mat","b_num");

figure; hold on; box off;
cmap_num = jet(N);
for j = 1:N
    plot(r_vec, b_num(:,j), 'LineWidth', 1.5, 'Color', cmap_num(j,:));
end
set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([0, 1]);
yticks([0 0.2 0.4 0.6 0.8 1]);
yticklabels({' ','0.2','0.4','0.6','0.8','1'});
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
ylabel('b(r)'); xlabel('r');
title('Numerical solution: Parameter \nu from 10^{-3} to 10^{2}');
subtitle(['\lambda = ', num2str(lambda), ',  \eta = ', num2str(eta)]);

P_dbl = 2*lambda/3;
Q_dbl = lambda*ones(N,1) - (((q_vec).^2)*eta)/2;
R_dbl = 1;

W_dbl = 3*eta*(q_vec).^2 .* P_dbl;
E_dbl = Q_dbl.^2 + W_dbl .* R_dbl;
F_dbl = (3/2).*(Q_dbl.*R_dbl - W_dbl./q_vec).*W_dbl + Q_dbl.^3;

pk = [0.5 + (sqrt(3)/2)*1i, 0.5 - (sqrt(3)/2)*1i, -1];

G_dbl = zeros(N, 3);
for m = 1:3
    G_dbl(:, m) = pk(m) .*(sqrt( (F_dbl - E_dbl.^(3/2)).*(F_dbl + E_dbl.^(3/2)) ) + F_dbl).^(1/3);
end

alpha_dbl = zeros(N, 3);
for m = 1:3
    alpha_dbl(:, m) = (1./W_dbl) .* (Q_dbl - G_dbl(:,m) - E_dbl./G_dbl(:,m));
end
alpha_opt = alpha_dbl(:, 1);

b_DblExpSer = zeros(N, N);

for j = 1:N
    q_j  = q_vec(j);

    if q_j < q_vec(60)
        al_j = alpha_opt(j);
        use_real_output = true;
    else
        al_j = alpha_opt(j);
        use_real_output = false;
    end

    exp_j = exp(-al_j .* r_vec);
    b_tmp = 1 + q_j .* (exp_j - 1)./r_vec ...
          - (q_j^3).*eta .* ((exp_j.*(1 + al_j.*r_vec) - 1)) ./ (r_vec.^2);

    if use_real_output
        b_DblExpSer(:, j) = real(b_tmp);
    else
        b_DblExpSer(:, j) = b_tmp;
    end
end

err_DblExpSer      = abs(b_DblExpSer - b_num);
mean_err_DblExpSer = sum(err_DblExpSer, 1)' ./ N;
log_err_DblExpSer  = log10(mean_err_DblExpSer);

figure; hold on; box on;
plot(nu_vec, log_err_DblExpSer, '--', 'Color', 'm', 'LineWidth', 2.5);
set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([-5.5, -0.5]);
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
yticks([-5.5 -5 -4 -3 -2 -1 -0.5]);
yticklabels({' ','-5','-4','-3','-2','-1',' '});
xlabel('\nu'); ylabel('log_{10}||b_{num} - b_{approx}||');
title('Analytical approximation: DblExpSer. Parameter \nu from 10^{-3} to 10^{2}');
legend('DblExpSer', 'Location', 'southeast');

S_expvar = zeros(N,1);
alpha_expvar = zeros(N,1);
ExpVar = zeros(N,N);

for i = 1:N
     q1 = q_vec(i);
     S_expvar(i,1) = (1+2*q1*eta)/3;
     alpha_expvar(i,1) = (sqrt(S_expvar(i,1).^2 + lambda./ q1) - S_expvar(i,1))./lambda;
    for j = 1:N
        ExpVar(j,i) = 1 + q1.*((exp(-alpha_expvar(i,1).*r_vec(j)) - 1)./r_vec(j));
    end
end

err_ExpVar     = abs(ExpVar - b_num);
mean_err_ExpVar = sum(err_ExpVar, 1)' ./ N;
log_err_ExpVar  = log10(mean_err_ExpVar);

figure; hold on; box on;
plot(nu_vec, log_err_ExpVar, '--', 'Color', 'b', 'LineWidth', 2.5);
set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([-5.5, -0.5]);
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
yticks([-5.5 -5 -4 -3 -2 -1 -0.5]);
yticklabels({' ','-5','-4','-3','-2','-1',' '});
xlabel('\nu'); ylabel('log_{10}||b_{num} - b_{approx}||');
title(['\lambda = ', num2str(lambda), ',  \eta = ', num2str(eta)]);
legend('ExpVar', 'Location', 'southeast');

S_expglobal =zeros(N,1);
alpha_expglobal =zeros(N,1);
ExpGlobal = zeros(N,N);

for i = 1:N
    q1 = q_vec(i);
    S_expglobal(i,1) = log(3/2) + q1*eta*log(4/3);
    alpha_expglobal(i,1) = (sqrt(S_expglobal(i,1)^2 + lambda/q1) - S_expglobal(i,1))/lambda;
    for j = 1:N
        ExpGlobal(j,i) = 1 + q1*((exp(-alpha_expglobal(i,1).*r_vec(j)) - 1)./r_vec(j));
    end
end

err_ExpGlobal   = abs(ExpGlobal - b_num);
mean_ExpGlobal  = sum(err_ExpGlobal, 1)' ./ N;
log_err_ExpGlobal   = log10(mean_ExpGlobal);

figure; hold on; box on;
plot(nu_vec, log_err_ExpGlobal, '-', 'Color', 'b', 'LineWidth', 2.5);
set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([-5.5, -0.5]);
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
yticks([-5.5 -5 -4 -3 -2 -1 -0.5]);
yticklabels({' ','-5','-4','-3','-2','-1',' '});
xlabel('\nu'); ylabel('log_{10}||b_{num} - b_{approx}||');
title(['\lambda = ', num2str(lambda), ',  \eta = ', num2str(eta)]);
legend('ExpGlobal', 'Location', 'southeast');

pk = [1/2 + (sqrt(3)/2)*1i, 1/2 - (sqrt(3)/2)*1i,-1];
P_Dbl_Exp_Var = zeros(N,1);
Q_Dbl_Exp_Var = zeros(N,1);
R_Dbl_Exp_Var = zeros(N,1);
W_Dbl_Exp_Var = zeros(N,1);
E_Dbl_Exp_Var = zeros(N,1);
F_Dbl_Exp_Var = zeros(N,1);

G_Dbl_Exp_Var = zeros(N,3);
alpha_Dbl_Exp_Var = zeros(N,3);

for i = 1:N
    q = q_vec(i);
    P_Dbl_Exp_Var(i,1) = lambda*(8*log(2) - 5) + ((4/3)*(q^2)*eta)*(1-eta*q)*(1-3*log(4/3));
    Q_Dbl_Exp_Var(i,1) = lambda + ((2/3)*(q^2)*eta)*(1-6*log(9/8) +(2*q*eta)*(1-6*log(4/3)));
    R_Dbl_Exp_Var(i,1) = (q*eta + 2)/3;
    W_Dbl_Exp_Var(i,1) = 3*eta*(q^2)*P_Dbl_Exp_Var(i,1);
    E_Dbl_Exp_Var(i,1) = (Q_Dbl_Exp_Var(i,1)^2) + W_Dbl_Exp_Var(i,1)*R_Dbl_Exp_Var(i,1);
    F_Dbl_Exp_Var(i,1) = (3/2)*(Q_Dbl_Exp_Var(i,1)*R_Dbl_Exp_Var(i,1) ...
        - (W_Dbl_Exp_Var(i,1)/q))*W_Dbl_Exp_Var(i,1)...
        + Q_Dbl_Exp_Var(i,1)^3;

end
    for j = 1:3
    G_Dbl_Exp_Var(:,j) = pk(j)*(sqrt((F_Dbl_Exp_Var).^2 - (E_Dbl_Exp_Var).^3) + F_Dbl_Exp_Var).^(1/3);
    alpha_Dbl_Exp_Var(:,j) = (1./W_Dbl_Exp_Var).*(Q_Dbl_Exp_Var - G_Dbl_Exp_Var(:,j) - (E_Dbl_Exp_Var./G_Dbl_Exp_Var(:,j)));
    end

ALPHA_Dbl_Exp_Var = alpha_Dbl_Exp_Var(:,1);
DblExpVar = zeros(N,N);

for k = 1:N
    q = q_vec(k);
    for i = 1:N
    DblExpVar(i,k) = 1 + q*((exp(-ALPHA_Dbl_Exp_Var(k).*r_vec(i)) - 1)./r_vec(i))-...
    (q^3)*eta*((exp(-ALPHA_Dbl_Exp_Var(k).*r_vec(i)).*(1+ALPHA_Dbl_Exp_Var(k).*r_vec(i)) - 1)./(r_vec(i).^2));
    end
end

err_DblExpVar   = abs(DblExpVar - b_num);
mean_DblExpVar  = sum(err_DblExpVar, 1)' ./ N;
log_err_DblExpVar   = log10(mean_DblExpVar);

figure; hold on; box on;
plot(nu_vec, log_err_DblExpVar, '--', 'Color', 'm', 'LineWidth', 2.5);
set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([-5.5, -0.5]);
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
yticks([-5.5 -5 -4 -3 -2 -1 -0.5]);
yticklabels({' ','-5','-4','-3','-2','-1',' '});
xlabel('\nu'); ylabel('log_{10}||b_{num} - b_{approx}||');
title(['\lambda = ', num2str(lambda), ',  \eta = ', num2str(eta)]);
legend('ExpGlobal', 'Location', 'southeast');

pk = [1/2 + (sqrt(3)/2)*1i, 1/2 - (sqrt(3)/2)*1i,-1];
P_Dbl_Exp_Global = zeros(N,1);
Q_Dbl_Exp_Global = zeros(N,1);
R_Dbl_Exp_Global = zeros(N,1);
W_Dbl_Exp_Global = zeros(N,1);
E_Dbl_Exp_Global = zeros(N,1);
F_Dbl_Exp_Global = zeros(N,1);

G_Dbl_Exp_Global = zeros(N,3);
alpha_Dbl_Exp_Global = zeros(N,3);

for i = 1:N
   q = q_vec(i);
   P_Dbl_Exp_Global(i,1) = 2*lambda*(1-log(2)) + q^2*eta*(1-q*eta)*(log(3) - 1);
   Q_Dbl_Exp_Global(i,1) = lambda - 2*q^2*eta*(1-log(81/32) + 2*q*eta*log(9/8));
   R_Dbl_Exp_Global(i,1) = q*eta + 2*(1-q*eta)*log(3/2);
   W_Dbl_Exp_Global(i,1) = 3*eta*q^2*P_Dbl_Exp_Global(i,1);
   E_Dbl_Exp_Global(i,1) = Q_Dbl_Exp_Global(i,1)^2 + W_Dbl_Exp_Global(i,1)*R_Dbl_Exp_Global(i,1);
   F_Dbl_Exp_Global(i,1) = (3/2)*(Q_Dbl_Exp_Global(i,1)*R_Dbl_Exp_Global(i,1) ...
    - W_Dbl_Exp_Global(i,1)/q)*W_Dbl_Exp_Global(i,1) + Q_Dbl_Exp_Global(i,1)^3;
end

for j = 1:3
    G_Dbl_Exp_Global(:,j) = pk(j)*(sqrt((F_Dbl_Exp_Global).^2 ...
        - (E_Dbl_Exp_Global).^3) ...
        + F_Dbl_Exp_Global).^(1/3);
    alpha_Dbl_Exp_Global(:,j) = (1./W_Dbl_Exp_Global).*(Q_Dbl_Exp_Global- G_Dbl_Exp_Global(:,j) - ( E_Dbl_Exp_Global./G_Dbl_Exp_Global(:,j) ) );
end
ALPHA_Dbl_Exp_Global =  alpha_Dbl_Exp_Global(:,1);

DblExpGlobal = zeros(N,N);

for k = 1:N
    q = q_vec(k);
    for i = 1:N
    DblExpGlobal(i,k) = 1 + q*((exp(-ALPHA_Dbl_Exp_Global(k)*r_vec(i)) - 1)/r_vec(i))-...
    (q^3)*eta*((exp(-ALPHA_Dbl_Exp_Global(k).*r_vec(i)).*(1+ALPHA_Dbl_Exp_Global(k).*r_vec(i)) - 1)./(r_vec(i).^2));
    end
end

err_DblExpGlobal   = abs(DblExpGlobal - b_num);
mean_DblExpGlobal  = sum(err_DblExpGlobal, 1)' ./ N;
log_err_DblExpGlobal   = log10(mean_DblExpGlobal);

figure; hold on; box on;
plot(nu_vec, log_err_DblExpGlobal, '-', 'Color', 'm', 'LineWidth', 2.5);
set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([-5.5, -0.5]);
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
yticks([-5.5 -5 -4 -3 -2 -1 -0.5]);
yticklabels({' ','-5','-4','-3','-2','-1',' '});
xlabel('\nu'); ylabel('log_{10}||b_{num} - b_{approx}||');
title(['\lambda = ', num2str(lambda), ',  \eta = ', num2str(eta)]);
legend('ExpGlobal', 'Location', 'southeast');

figure; hold on; box on;
plot(nu_vec, log_err_DblExpSer, '--', 'Color',  'm', 'LineWidth', 2.5);
plot(nu_vec, log_err_ExpVar,       '--', 'Color',  'b', 'LineWidth', 2.5);
plot(nu_vec, log_err_ExpGlobal,    '-', 'Color',  'b', 'LineWidth', 2.5);
plot(nu_vec, log_err_DblExpVar, ':.', 'Color', 'm', 'LineWidth', 2.5);
plot(nu_vec, log_err_DblExpGlobal,    '-',  'Color',  'm', 'LineWidth', 1.5);

set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([-5.5, -0.5]);
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
yticks([-5.5 -5 -4 -3 -2 -1 -0.5]);
yticklabels({' ','-5','-4','-3','-2','-1',' '});
xlabel('\nu'); ylabel('log_{10}||b_{num} - b_{approx}||');
title({'Error comparison with DblExpSer, ExpVar, ExpGlobal','DblExpVar and DblExpGlobal'});
subtitle(['\lambda = ', num2str(lambda), ',  \eta = ', num2str(eta)]);
legend('DblExp-Ser','ExpVar','ExpGlobal','DblExp-Var','DblExp-Global','Location', 'southeast');
