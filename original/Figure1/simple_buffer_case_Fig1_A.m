%% 
clear, close all; clc;
%% 参数
R_max = 1e+2;
r_eps = 1e-15;         
N = 5e+4;
b0 = 0.131;
lambda = 0.1;
nu = 0.1;
eta = 1;
q = 1/(eta + nu);
r = linspace(r_eps, R_max, N)';
h = r(2) - r(1);
P_dbl = (2*lambda) /3;
Q_dbl = lambda - (q^2*eta) / 2;
R_dbl = 1;

% Pade2 parameter
p = q*nu;
R = 24 * lambda - 3 * p * q * (p - 2);
K = 8 * lambda^2 - 2 * lambda * q * (2*p^2 - 5*p - 2) - p * q^2 * (p - 2);
Q = K + 24 * lambda^2 - 4 * lambda * p * q * (2*p - 3);
V = 2 * lambda * (p + 6) + p * q * (p - 3) * (p - 2);
F = 6 * lambda^2 - lambda * q * (p + 3) * (p - 2) - eta * p * q^3 * (2*p - 3);
J = 24 * lambda^2 - 2 * lambda * p * q * (5*p - 8) - eta * p^2 * q^3 * (p - 2);
E = 12 * lambda^2 + eta * p * q^2 * (8 * lambda + eta * p * q^2);
H = p^2 * q^4 * (p-2)^2 + 1600 * lambda^4 ...
    - 104 * lambda^3 * q * (p-2) * (9*p + 4) ...
    + 4 * lambda^2 * q^2 * (39*p^4 - 87*p^3 - 63*p^2 + 184*p + 4) ...
    - 4 * lambda * p * q^3 * (p-2) * (6*p^3 - 21*p^2 + 20*p + 2);
W = 12000 * lambda^5 ...
    - 4 * lambda^4 * q * (2687*p^2 - 2966*p - 553) ...
    + 4 * lambda^3 * q^2 * (729*p^4 - 1626*p^3 + 486*p^2 + 421*p + 70) ...
    - lambda^2 * q^3 * (243*p^6 - 810*p^5 + 609*p^4 + 210*p^3 - 21*p^2 - 240*p - 4) ...
    - 2 * lambda * eta * p * q^5 * (p-2) * (20*p^2 - 13*p - 1) ...
    + p^2 * eta^2 * q^7;
% G 中 sqrt(-3 * q * W) 注意定义域，这里直接保留表达式
G = Q^3 + 2 * lambda * R * (3 * Q * F - 9 * lambda * E * R + V * sqrt(-3 * q * W));

B2 = (q / R) * (Q + G^(1/3) + H / G^(1/3));
B1 = (1 / V) * ((B2 * (B2 * R - 3 * q * K)) / (6 * lambda * q^2) - J);
A1 = B1 - q;
A2 = B2 - q * (B1 - eta * q^2);

%EXP-PADÉ Parameter APPENDIX B
P_ExpPade = 6 + (eta * q^2 * (1 - eta * q)) / lambda;
Q_ExpPade = (1 / P_ExpPade) * (6 + (eta * q^2) / (2 * lambda));
R_ExpPade = Q_ExpPade * (1 - (2 * lambda) / q) + (2 * eta^2 * q^2) / P_ExpPade;
U_ExpPade = R_ExpPade / 3 - Q_ExpPade^2 / 2;
V_ExpPade = Q_ExpPade * (R_ExpPade - Q_ExpPade^2) + (12 * lambda) / (q * P_ExpPade);
W_ExpPade = Q_ExpPade^2 * (R_ExpPade / 3 - Q_ExpPade^2 / 4) ...
            + 4 * lambda * (6 * (Q_ExpPade + lambda / q) - eta^2 * q^2) / (3 * q * P_ExpPade);
E_ExpPade = V_ExpPade^2 * (V_ExpPade^2 / 4 + U_ExpPade * (U_ExpPade^2 - 3 * W_ExpPade)) ...
            - W_ExpPade * (W_ExpPade - 3 * U_ExpPade^2)^2;
G_ExpPade = V_ExpPade^2 / 2 + U_ExpPade * (U_ExpPade^2 - 3 * W_ExpPade) + E_ExpPade^(1/2);
H_ExpPade = (W_ExpPade + U_ExpPade^2) / G_ExpPade^(1/3) +  G_ExpPade^(1/3) - 2 * U_ExpPade;

alpha_ExpPade = (1 / (2 * lambda)) * (H_ExpPade^(1/2) - Q_ExpPade + (2 * V_ExpPade * H_ExpPade^(-1/2) - H_ExpPade - 6 * U_ExpPade)^(1/2));
beta_ExpPade = eta * q^3 / (q * alpha_ExpPade * (lambda * alpha_ExpPade + 1) - 1);

%% 初始猜测：分段构造
r_switch = 3.00;            % 切换点（可调）
b_init = zeros(N,1);
for i = 1:N
    ri = r(i);
    if ri < r_switch
        % 小 r：Frobenius 展开
        b_init(i) = b0 + (b0/(2*lambda))*ri + ...
            (((b0-1)*(nu*b0+eta) + b0/(2*lambda))/(6*lambda))*ri^2;
    else
        % 大 r：渐近展开 
        x = 1/ri;
        b_init(i) = 1 - q*x + eta*q^3*x^2 + eta*(1-2*q*eta)*q^4*x^3 + ...
            eta*(2*lambda + 5*eta*q^2*(eta*q-1))*q^4*x^4;
    end
end
b = b_init;
% % 可选：在 r_switch 附近做平滑连接（例如线性加权）
% if r_switch > r_min && r_switch < R_max
%     idx = find(r >= r_switch, 1, 'first');
%     if idx > 1 && idx <= N
%         w = (r(idx) - r_switch) / (r(idx) - r(idx-1));  % 实际上 r_switch 落在两节点间
%         % 这里简化：直接取两节点的平均值，或不做平滑（作为初始猜测足够）
%     end
% end

% 牛顿迭代参数
tol = 1e-8;
max_iter = 25;

for iter = 1:max_iter
    [F, J] = compute_FJ_Fig1(b, r, h, N, lambda, nu, eta, b0);
    b = b - J \ F;
    res = norm(F, inf);
    fprintf('Iter %d: residual = %.3e\n', iter, res);
    if res < tol
        fprintf('Converged.\n');
        break;
    end
end
%% RBA2
% b_RPA 关于b的一元二次方程有两个根，大一点的是这个，带进去的话算出来是另一种不知意义的解
%b_RPA = ((nu - eta - 1./r) + sqrt((eta - nu + 1./r).^2 + 4*nu*(eta + 1./r)))/(2*nu);
% 最后将b_RPA换成了RBA
b_RPA = 1 - 1/(2*q*nu)*(1+q./r - sqrt( (1+q./r).^2 - 4*nu*q^2./r ));
RBA2 = b_RPA + 2*lambda*eta*( (1+r/q).^2 - 4*nu.*r ).^(-2);

%% Pade
Frac_under = q ./ (r + (sqrt(q*(q+8*lambda)) + q)/2);
Pade = 1 - Frac_under;

%% Pade2
Pade2 = (r.^2 + A1.*r + A2) ./ (r.^2 + B1.*r + B2);

%% Exp-Ser
S = 1/2;
alpha_expser = (sqrt(S^2 + lambda/q) - S) / lambda;
expser_funcform_3 = 1 + q*((exp(-alpha_expser*r)-1)./r);

%% EXP-PADÉ
Exp_Pade = 1 + q*(exp(-alpha_ExpPade*r)-1)./r + q^3*eta*1./(beta_ExpPade + r.^2);

%% DbleExpSer
syms alpha
DbleExpSer = eta*q^2*P_dbl*alpha^3 - Q_dbl*alpha^2 - R_dbl*alpha + 1/q;
DbleExpSer_coeff = [eta*q^2*P_dbl,-Q_dbl,-R_dbl,1/q];
alpha_roots = roots(DbleExpSer_coeff);
new_alpha = 1.2558+0.926i; %拿到alpha之后不取实部
b_DbleExpSer0 = real(1-q*new_alpha + 1/2*q^3*new_alpha^2);

b_funcform_2 = 1 + (q.*(exp(-new_alpha*r)-1))./r ...% 用functional form计算完之后再取实部
    - q^3*eta.*((exp(-new_alpha*r).*(1+new_alpha*r)-1)./r.^2); 
%% 绘图
figure;
hold on; box off;
plot(r, b, 'Color', [0.75 0.75 0.75], 'LineWidth', 2);
plot(r, RBA2, '--', 'Color','r', 'LineWidth',2);
plot(r, Pade, '-', 'Color',[46/256,139/256,87/256], 'LineWidth',2);
plot(r, real(Pade2), '--', 'Color', [46/256,139/256,87/256], 'LineWidth',2);
plot(r, expser_funcform_3, '-', 'Color','k', 'LineWidth',2);
plot(r, real(Exp_Pade), '--', 'Color','k', 'LineWidth',2);
plot(r, real(b_funcform_2), '--', 'Color','m', 'LineWidth',2);

ylabel('b(r)');
xlim([1e-3, 1e+2]);
xticks([1e-3 1e-2 1e-1 1 10 100]);
ylim([0.089 1]);
yticks([0 0.2 0.4 0.6 0.8 1]);
yticklabels({' ','0.2','0.4','0.6','0.8','1'});
xlabel('r (distance)');
set(gca, 'XScale', 'log');
title('Numerical solution and Analytical approximations');
subtitle('\lambda = 0.1,\nu = 0.1, \eta = 1')
legend('Numerical', 'RBA2','Padé','Padé2','Exp-Ser','EXP-Padé','DblExp-Ser','Location', 'southeast');
exportgraphics(gcf, 'simple_buffer_case_Fig1_A.png', ...
    'Resolution', 300, 'BackgroundColor', 'white');