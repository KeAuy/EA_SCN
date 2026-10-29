clear; close all; clc;

%% System paramter
lambda   = 0.02;   
eta      = 1.0;      
N        = 100;    

tol      = 1e-7;   
max_iter = 50;     
nu_vec = logspace(-3, 2, N)'; 
%% Intermediate 

k_idx  = (1:N)';
r_vec  = 10.^(-3 + 5.*k_idx./N);   % N×1
q_vec  = 1 ./ (eta + nu_vec);      % N×1
x_vec  = 1 ./ r_vec;               % N×1

% 非线性网格差分的 left-step 和 right step
alpha_r = 10^(5/N);                 % r_{k+1}/r_k = alpha_r
h_L = r_vec .* (1 - 1/alpha_r);    % left-step  h^- = r_k - r_{k-1},  N×1
h_R = r_vec .* (alpha_r - 1);      % right-step h^+ = r_{k+1} - r_k,  N×1

%% shooting method
temp  = load("b0_Fig4_A.mat");
b_r0  = temp.b0;   
%% Asymptotic expansion 

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

%% switch rad 
% 选取使两种展开相对误差最小的 r 作为拼接点：
rel_disc   = abs(b_frob - b_asymp) ./ (abs(b_frob) + abs(b_asymp));  % N×N
[~, idx_switch] = min(rel_disc, [], 1);   % 1×N 
r_switch   = r_vec(idx_switch);           % 1×N 切换半径

%%  b_init 
% r < r_switch 时用 Frobenius，否则用无穷远处展开
b_init = zeros(N, N);
 
for j = 1:N
    mask = r_vec < r_switch(j);
    %写成 Characteristic function 的形式
    b_init(:, j) = double(mask).*b_frob(:, j) ...       
                 + double(~mask).*b_asymp(:, j);
end
%% Newton Iteration 
b_num       = b_init;          % N×N
res_vec     = zeros(N, 1);     
failed_cols = zeros(2, N);     % 未收敛记录：[列索引; ν值]

for j = 1:N
    b0j        = b_r0(j);
    b_col_init = b_init(:, j);   
    converged  = false;

    % Step 1：Newton  
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

    % Step 2： 引入松弛项
    if ~converged
        fprintf('  Newton False，j=%d，引入松弛项。\n', j);
        b_num(:, j) = b_col_init;   

        for iter = 1:max_iter
            F_col   = compute_F(b0j, b_num(:,j), r_vec, h_L, h_R, N, lambda, nu_vec(j), eta);
            res_old = norm(F_col, inf);
            if res_old < tol
                fprintf(' 松弛法迭代结果收敛，j=%d\n', j);
                converged = true;
                break;
            end
            J_col = compute_J(b_num(:,j), r_vec, h_L, h_R, N, lambda, nu_vec(j), eta);
            delta = -J_col \ F_col;

            % 引入松弛项
            omega_grid = 0.1:0.05:1.0;
            best_res   = inf;
            best_b     = b_num(:, j);
            for omega = omega_grid
                b_trial   = b_num(:, j) + omega * delta;
                F_trial   = compute_F(b0j, b_trial, r_vec, h_L, h_R, N, lambda, nu_vec(j), eta);
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
save("b_num_Fig4_A.mat","b_num");   
%%  解析近似：RBA2
b_RBA  = zeros(N, N);
b_RBA2 = zeros(N, N);

for j = 1:N
    nu_j = nu_vec(j);
    q_j  = q_vec(j);
    for k = 1:N
        r_k = r_vec(k);
        disc_rba     = (1 + q_j/r_k)^2 - 4*nu_j*(q_j^2)/r_k;
        b_RBA(k, j)  = 1 - (1/(2*nu_j*q_j)) * (1 + q_j/r_k - sqrt(disc_rba));
        b_RBA2(k, j) = b_RBA(k, j) + 2*lambda*eta / ((1 + r_k/q_j)^2 - 4*nu_j*r_k)^2;
    end
end
%% Error: RBA2
err_RBA2       = abs(b_RBA2 - b_num);                 
mean_err_RBA2  = sum(err_RBA2, 1)' ./ N;             
log_err_RBA2   = log10(mean_err_RBA2);

figure; hold on; box on;
plot(nu_vec, log_err_RBA2, '--', 'Color', 'r', 'LineWidth', 2.5);
set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([-5.5, -0.5]);
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
yticks([-5.5 -5 -4 -3 -2 -1 -0.5]);
yticklabels({' ','-5','-4','-3','-2','-1',' '});
xlabel('\nu');ylabel('log_{10}||b_{num} - b_{approx}||');
title('RBA2 Error: \lambda = 0.02,\eta = 1');
legend('RBA2', 'Location', 'southeast');

%% 解析近似: Padé 2 
b_Pade2 = zeros(N,N);
p_vec  = q_vec .* nu_vec;  

Pc_R = (24*lambda) - 3.*p_vec.*q_vec.*(p_vec - 2);
Pc_K = 8*(lambda^2) ...
       - 2*lambda.*q_vec.*(2*(p_vec).^2 - 5*p_vec - 2) ...
       - p_vec.*((q_vec).^2).*(p_vec - 2);
Pc_Q = Pc_K + 24*(lambda^2) - 4*lambda.*p_vec.*q_vec.*(2*p_vec - 3);
Pc_V = 2*lambda.*(p_vec + 6) + p_vec.*q_vec.*(p_vec - 3).*(p_vec - 2);
Pc_F = 6*(lambda^2) ...
       - lambda.*q_vec.*(p_vec + 3).*(p_vec - 2) ...
       - eta.*p_vec.*(q_vec.^3).*(2*p_vec - 3);
Pc_J = 24*(lambda^2) ...
       - 2*lambda.*p_vec.*q_vec.*(5*p_vec - 8) ...
       - eta.*(p_vec.^2).*(q_vec.^3).*(p_vec - 2);
Pc_E = 12*(lambda^2) + eta.*p_vec.*(q_vec.^2).*(8*lambda + eta.*p_vec.*(q_vec.^2));
Pc_H = (p_vec.^2).*(q_vec.^4).*((p_vec - 2).^2) + 1600*(lambda^4) ...
       - 104*(lambda^3).*q_vec.*(p_vec - 2).*(9*p_vec + 4) ...
       + 4*(lambda^2).*(q_vec.^2) ...
         .*(39*(p_vec.^4) - 87*(p_vec.^3) - 63*(p_vec.^2) + 184*p_vec + 4) ...
       - 4*lambda.*p_vec.*(q_vec.^3).*(p_vec - 2) ...
         .*(6*(p_vec.^3) - 21*(p_vec.^2) + 20*p_vec + 2);
Pc_W = 12000*(lambda^5) ...
       - 4*(lambda^4).*q_vec.*(2687*(p_vec.^2) - 2966*p_vec - 553) ...
       + 4*(lambda^3).*(q_vec.^2) ...
         .*(729*(p_vec.^4) - 1626*(p_vec.^3) + 486*(p_vec.^2) + 421*p_vec + 70) ...
       - (lambda^2).*(q_vec.^3) ...
         .*(243*(p_vec.^6) - 810*(p_vec.^5) + 609*(p_vec.^4) ...
            + 210*(p_vec.^3) - 21*(p_vec.^2) - 240*p_vec - 4) ...
       - 2*lambda*eta.*p_vec.*(q_vec.^5).*(p_vec - 2).*(20*(p_vec.^2) - 13*p_vec - 1) ...
       + (p_vec.^2).*(eta^2).*(q_vec.^7);

r_qw = abs(-3.*q_vec .* Pc_W);
theta_qw = angle(-3.*q_vec .* Pc_W);
sqrt_qw = exp((log(r_qw) + 1i*theta_qw)/2);

Pc_G   = Pc_Q.^3 + 2*lambda.*Pc_R ...
          .* (3*Pc_Q.*Pc_F - 9*lambda.*Pc_E.*Pc_R + Pc_V.*(sqrt_qw));

B2 = real((q_vec ./ Pc_R) .* (Pc_Q + (Pc_G).^(1/3) + (Pc_H ./ (Pc_G).^(1/3) )));
B1 = (1./Pc_V) .* (( (B2.*(B2.*Pc_R - 3.*q_vec.*Pc_K)) ./ (6*lambda*(q_vec).^2) )- Pc_J );
A1 = B1 - q_vec;
A2 = B2 - q_vec .* (B1 - eta*(q_vec.^2));

% ============================================================
% 逐列求解代数方程组 (η = 1 情形)
%   2λB₁² + B₁B₂ - B₂²/q - 2λq²B₁ - (q²+2λ)B₂ = 0               ... (1)
%   (1+q)B₁² - 2B₁B₂/q + 2(3λ-q³)B₁ + (1-q)B₂ = q²(6λ+q²(1-q))  ... (2)
% ============================================================

B1_vec = zeros(N,1);
B2_vec = zeros(N,1);

options = optimoptions('fsolve','Display','none','TolFun',1e-14,...
    'TolX',1e-14,'MaxIterations',500);
flag = zeros(N,1);

for i = 1:N
    q = q_vec(i);

    parm_b_pade = @(B) [2*lambda*B(1)^2 + B(1)*B(2)-B(2)^2/q ...
        - 2*lambda*q^2*B(1) - (q^2 + 2*lambda)*B(2);
        (1+q)*B(1)^2 - 2*B(1)*B(2)/q + 2*(3*lambda-q^3)*B(1) ...
        + (1-q)*B(2) - q^2*(6*lambda + q^2*(1-q))];
    
    B0 = [B1(i); B2(i)];
    [B, ~, exitflag] = fsolve(parm_b_pade, B0, options);
    if exitflag <= 0 && i > 1
        B0 = [B1(i-1); B2(i-1)];
        [B,~,~] = fsolve(parm_b_pade, B0, options);
    end
    B1_vec(i) = B(1);
    B2_vec(i) = B(2);
end


A1_vec = B1_vec - q_vec;
A2_vec = B2_vec - q_vec.*(B1_vec -eta* q_vec.^2);

for j=1:N
    b_Pade2(:,j) = (r_vec.^2 + A1_vec(j).*r_vec + A2_vec(j)) ...
    ./ (r_vec.^2 + B1_vec(j).*r_vec + B2_vec(j));
end
%% Plot: Pade2
figure; hold on; box off;
cmap_pade = jet(N);
for j = 1:N
    plot(r_vec, b_Pade2(:,j), 'LineWidth', 1.5, 'Color', cmap_pade(j,:));
end
set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([0, 1]);
yticks([0 0.2 0.4 0.6 0.8 1]); yticklabels({' ','0.2','0.4','0.6','0.8','1'});
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
ylabel('b(r)'); xlabel('r');
title('Analytical approximation: Padé2,  \nu \in [10^{-3}, 10^2]');
subtitle(['\lambda = ', num2str(lambda), ',  \eta = ', num2str(eta)]);
hold off;
%% Error:Pade2
err_Pade2      = abs(b_Pade2 - b_num);
mean_err_Pade2 = sum(err_Pade2, 1)' ./ N;
log_err_Pade2  = log10(mean_err_Pade2);

figure; hold on; box on;
plot(nu_vec, log_err_Pade2, '--', 'Color', 'g', 'LineWidth', 2.5);
set(gca, 'XScale', 'log');
xlim([1e-3, 1e2]); ylim([-5.5, -0.5]);
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
yticks([-5.5 -5 -4 -3 -2 -1 -0.5]);
yticklabels({' ','-5','-4','-3','-2','-1',' '});
xlabel('\nu'); ylabel('log_{10}||b_{num} - b_{approx}||');
title(['\lambda = ', num2str(lambda), ',  \eta = ', num2str(eta)]);
legend('Padé2', 'Location', 'southeast');
%% 解析近似: DblExpSer

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
%%
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
%% Error: DblExpSer
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

%% 解析近似: Exp-Var
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
%% Error:Exp-Var
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
%% 解析近似：Exp-Global Parameter
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
%% Error:Exp-Global

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


%%  解析近似：DblExp-var
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
%% Error: DblExp-var
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

%% DblExp-Global 
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

%% Error: DblExp-Global
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


%% 误差汇总对比
figure; hold on; box on;
plot(nu_vec, log_err_RBA2,      '--', 'Color',  'r', 'LineWidth', 2.5);
plot(nu_vec, log_err_Pade2,     '--', 'Color',  'g', 'LineWidth', 2.5);
plot(nu_vec, log_err_DblExpSer, '--', 'Color',  'm', 'LineWidth', 2.5);
plot(nu_vec, log_err_ExpVar,       '--', 'Color',  'b', 'LineWidth', 2.5);
plot(nu_vec, log_err_ExpGlobal,    '-', 'Color',  'b', 'LineWidth', 2.5);
plot(nu_vec, log_err_DblExpVar, ':.', 'Color', 'm', 'LineWidth', 2.5);
plot(nu_vec, log_err_DblExpGlobal,    '-',  'Color',  'm', 'LineWidth', 1.5);

set(gca, 'XScale', 'log');

% 实心点
[log_err_DblExpSer_min, log_err_DblExpSer_min_idx] = min(log_err_DblExpSer);   
log_err_DblExpSer_min_xpoint = nu_vec(log_err_DblExpSer_min_idx);
plot(log_err_DblExpSer_min_xpoint, log_err_DblExpSer_min, 'o', 'MarkerSize', 8, ...
     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'm');  

xlim([1e-3, 1e2]); ylim([-5.5, -0.5]);
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2]);
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
yticks([-5.5 -5 -4 -3 -2 -1 -0.5]);
yticklabels({' ','-5','-4','-3','-2','-1',' '});
xlabel('\nu'); ylabel('log_{10}||b_{num} - b_{approx}||');
title({'Error comparison with RBA2, Padé2, DblExp-Ser, ExpVar and', ...
       'ExpGlobal, DblExp-Var, DblExp-Global'});
subtitle(['\lambda = ', num2str(lambda), ',  \eta = ', num2str(eta)]);
legend('RBA2', 'Padé2','DblExp-Ser','ExpVar','ExpGlobal','DblExp-Var','DblExp-Global','Location', 'southeast');