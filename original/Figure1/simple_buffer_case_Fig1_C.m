clear,close all;clc;

R_max = 1e+2;
r_eps = 1e-9;
N = 5e+4;

lambda = 1;
nu = 10;
eta = 1;
q = 1/(eta + nu);
b0 = 0.7358;

r = linspace(r_eps,R_max, N)';
h = r(2) - r(1);

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
G = Q^3 + 2 * lambda * R * (3 * Q * F - 9 * lambda * E * R + V * sqrt(-3 * q * W));

B2 = (q / R) * (Q + G^(1/3) + H / G^(1/3));
B1 = (1 / V) * ((B2 * (B2 * R - 3 * q * K)) / (6 * lambda * q^2) - J);
A1 = B1 - q;
A2 = B2 - q * (B1 - eta * q^2);

S_ExpSer = 1/2;
alpha1C_ExpSer = (sqrt(S_ExpSer^2 + lambda/q) - S_ExpSer)/lambda;

P_DblExpSer = 2*lambda/3;
Q_DblExpSer = lambda - (q^2*eta)/2;
R_DblExpSer = 1;
W_DblExpSer = 3*eta*q^2*P_DblExpSer;
E_DblExpSer = Q_DblExpSer^2 + W_DblExpSer*R_DblExpSer;
F_DblExpSer = 3/2*(Q_DblExpSer*R_DblExpSer - W_DblExpSer/q)*W_DblExpSer + Q_DblExpSer^3;
pk = [1/2+sqrt(3)/2*i,1/2-sqrt(3)/2*i,-1];
G_k_DblExpSer = [];
alpha_DblExpSer = [];
for i = 1:3
     G_k_DblExpSer(i) = pk(i)*(sqrt(F_DblExpSer^2 - E_DblExpSer^3)+ F_DblExpSer)^(1/3);
end
for j = 1:3
     alpha_DblExpSer(j) = (1/W_DblExpSer)*(Q_DblExpSer - G_k_DblExpSer(j) - E_DblExpSer/G_k_DblExpSer(j));
end
new_alpha_DblExpSer = real(alpha_DblExpSer(1));

r_switch = 1.00;
b_init = zeros(N,1);
for i = 1:N
    ri = r(i);
    if ri < r_switch

        b_init(i) = b0 + (b0/(2*lambda))*ri + ...
            (((b0-1)*(nu*b0+eta) + b0/(2*lambda))/(6*lambda))*ri^2;
    else

        x = 1/ri;
        b_init(i) = 1 - q*x + eta*q^3*x^2 + eta*(1-2*q*eta)*q^4*x^3 + ...
            eta*(2*lambda + 5*eta*q^2*(eta*q-1))*q^4*x^4;
    end
end
b = b_init;

tol = 1e-10;
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

LIN_Fig1c = 1 + (q*((exp(-r/(sqrt(q*lambda))) - 1))./r);

Frac_under_1C = q ./ (r + (sqrt(q*(q+8*lambda)) + q)/2);
Pade_1C = 1 - Frac_under_1C;

Pade2 = (r.^2 + A1.*r + A2) ./ (r.^2 + B1.*r + B2);

ExpSer_1C = 1 +  q*(exp(-alpha1C_ExpSer*r)-1)./r;

DblExpSer_1C = 1 + (q.*(exp(-new_alpha_DblExpSer*r)-1))./r...
    - q^3*eta*((exp(-new_alpha_DblExpSer*r).*(1+new_alpha_DblExpSer*r)-1)./r.^2);

box off;
hold on;
set(gca,'XScale','log');xlim ([1e-3,1e+2]);
plot(r, b, 'Color',[0.75, 0.75, 0.75], 'LineWidth',2.5);
plot(r,LIN_Fig1c,'-.','Color','k','LineWidth',1.25);
plot(r, Pade_1C, '-', 'Color',[46/256,139/256,87/256], 'LineWidth',2.5);
plot(r, real(Pade2), '--', 'Color', [46/256,139/256,87/256], 'LineWidth',2.5);
plot(r, ExpSer_1C, '-', 'Color','k', 'LineWidth',2.5);
plot(r, DblExpSer_1C,'--', 'Color','m', 'LineWidth',2.5);
ylim([0.73,1]);
ax = gca;
labels = ax.YTickLabel;
labels{1} = '';
ax.YTickLabel = labels;
xlabel('r'); ylabel('b(r)');
xlabel('r(distance)')
yticks([0.73 0.75 0.8 0.85 0.9 0.95 1]);
yticklabels({' ','0.75','0.8','0.85','0.9','0.95','1'});
title('Numerical Solution and Analytic approximation');
subtitle('\lambda = 1, \nu = 10, \eta = 1')
xlabel('r(distance)')
legend('Numerical','LIN','Padé','Padé2','Exp-Ser','DblExp-Ser','Location', 'southeast');
exportgraphics(gcf, 'simple_buffer_case_Fig1_C.png', ...
    'Resolution', 300, 'BackgroundColor', 'white');
