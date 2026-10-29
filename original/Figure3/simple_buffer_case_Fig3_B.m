clear, close all; clc;

R_max = 1e+2;
r_eps = 1e-3;
N = 5e+4;
b0 = 0.369;

lambda = 0.1;
nu = 10;
eta = 1;

q = 1/(eta + nu);
r = linspace(r_eps, R_max, N)';
h = r(2) - r(1);
x = 1./r;

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

B2 = real((q / R) * (Q + G^(1/3) + H / G^(1/3)));
B1 = real((1 / V) * ((B2 * (B2 * R - 3 * q * K)) / (6 * lambda * q^2) - J));
A1 = B1 - q;
A2 = B2 - q * (B1 - eta * q^2);

S_expvar = (1+2*q*eta)/3;
alpha_expvar = (sqrt(S_expvar^2 + lambda/q) - S_expvar)/lambda;

S_expglobal = log(3/2) + q*eta*log(4/3);
alpha_expglobal = (sqrt(S_expglobal^2 + lambda/q) - S_expglobal)/lambda;

P_Dbl_Exp_Var = lambda*(8*log(2) - 5) + ((4/3)*(q^2)*eta)*(1-eta*q)*(1-3*log(4/3));
Q_Dbl_Exp_Var = lambda + ((2/3)*(q^2)*eta)*(1-6*log(9/8) +(2*q*eta)*(1-6*log(4/3)));
R_Dbl_Exp_Var = (q*eta + 2)/3;
W_Dbl_Exp_Var = 3*eta*(q^2)*P_Dbl_Exp_Var;
E_Dbl_Exp_Var = (Q_Dbl_Exp_Var^2) + W_Dbl_Exp_Var*R_Dbl_Exp_Var;
F_Dbl_Exp_Var = (3/2)*(Q_Dbl_Exp_Var*R_Dbl_Exp_Var ...
    - (W_Dbl_Exp_Var/q))*W_Dbl_Exp_Var...
    + Q_Dbl_Exp_Var^3;
pk = [1/2 + (sqrt(3)/2)*1i, 1/2 - (sqrt(3)/2)*1i,-1];
for i = 1:3
    G_Dbl_Exp_Var(i) = pk(i)*(sqrt((F_Dbl_Exp_Var)^2 ...
        - (E_Dbl_Exp_Var)^3) ...
        + F_Dbl_Exp_Var)^(1/3);
    alpha_Dbl_Exp_Var(i) = (1/W_Dbl_Exp_Var)...
        *(Q_Dbl_Exp_Var ...
        - G_Dbl_Exp_Var(i) ...
        - ( E_Dbl_Exp_Var/G_Dbl_Exp_Var(i) ) );
end
ALPHA_Dbl_Exp_Var = alpha_Dbl_Exp_Var(1);

P_Dbl_Exp_Global = 2*lambda*(1-log(2)) + q^2*eta*(1-q*eta)*(log(3) - 1);
Q_Dbl_Exp_Global = lambda - 2*q^2*eta*(1-log(81/32) + 2*q*eta*log(9/8));
R_Dbl_Exp_Global = q*eta + 2*(1-q*eta)*log(3/2);
W_Dbl_Exp_Global = 3*eta*q^2*P_Dbl_Exp_Global;
E_Dbl_Exp_Global = Q_Dbl_Exp_Global^2 + W_Dbl_Exp_Global*R_Dbl_Exp_Global;
F_Dbl_Exp_Global = (3/2)*(Q_Dbl_Exp_Global*R_Dbl_Exp_Global ...
    - W_Dbl_Exp_Global/q)*W_Dbl_Exp_Global + Q_Dbl_Exp_Global^3;

for i = 1:3
    G_Dbl_Exp_Global(i) = pk(i)*(sqrt((F_Dbl_Exp_Global)^2 ...
        - (E_Dbl_Exp_Global)^3) ...
        + F_Dbl_Exp_Global)^(1/3);
    alpha_Dbl_Exp_Global(i) = (1/W_Dbl_Exp_Global)...
        *(Q_Dbl_Exp_Global ...
        - G_Dbl_Exp_Global(i) ...
        - ( E_Dbl_Exp_Global/G_Dbl_Exp_Global(i) ) );
end
ALPHA_Dbl_Exp_Global =  alpha_Dbl_Exp_Global(1);

b_near0   = b0 ...
          + (b0/(2*lambda)) .* r ...
          + (((b0-1).*(nu'.*b0 + eta) + b0/(2*lambda)) ./ (6*lambda)) .* r.^2;

b_nearinf = 1 - q'.*x ...
          + eta.*(q.^3)'.*(x.^2) ...
          + eta.*(1 - 2*eta.*q').*(q.^4)'.*(x.^3) ...
          + eta.*((2*lambda) + 5*eta.*(q.^2)'.*(eta.*q' - 1)).*(q.^4)'.*(x.^4);

rel_diff     = abs(b_near0 - b_nearinf) ./ (abs(b_near0) + abs(b_nearinf) + 1e-30);
[~, idx_best] = min(rel_diff, [], 1);
r_switch      = r(idx_best);
[~, idx] = ismember(r_switch, r);

r_switch = r(idx);
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

tol = 1e-8;
max_iter = 25;

for iter = 1:max_iter
    [F, J] = compute_FJ_Fig3(b, r, h, N, lambda, nu, eta, b0);
    b = b - J \ F;
    res = norm(F, inf);
    fprintf('Iter %d: residual = %.3e\n', iter, res);
    if res < tol
        fprintf('Converged.\n');
        break;
    end
end

Pade2 = (r.^2 + A1.*r + A2) ./ (r.^2 + B1.*r + B2);

ExpVar = 1 + q*( (exp(-alpha_expvar.*r) - 1)./r );

ExpGlobal = 1 + q*( (exp(-alpha_expglobal.*r) - 1)./r );

DblExpVar = 1 + q*((exp(-ALPHA_Dbl_Exp_Var.*r) - 1)./r)-...
    (q^3)*eta*((exp(-ALPHA_Dbl_Exp_Var.*r).*(1+ALPHA_Dbl_Exp_Var.*r) - 1)./(r.^2));

DblExpGlobal = 1 + q*((exp(-ALPHA_Dbl_Exp_Global.*r) - 1)./r)-...
    (q^3)*eta*((exp(-ALPHA_Dbl_Exp_Global.*r).*(1+ALPHA_Dbl_Exp_Global.*r) - 1)./(r.^2));

figure;
hold on; box off;
plot(r, b, 'Color', [0.75 0.75 0.75], 'LineWidth', 2);
plot(r, real(Pade2), '--', 'Color', [46/256,139/256,87/256], 'LineWidth',2.5);
plot(r,ExpVar,'--', 'Color', [0,54/256,253/256], 'LineWidth',2);
plot(r,ExpGlobal,'-','Color',[0,54/256,253/256], 'LineWidth',2);
plot(r,real(DblExpVar),'--','Color','m','LineWidth',2);
plot(r,real(DblExpGlobal),'-','Color','m','LineWidth',2);
ylabel('b(r)');
xlim([1e-3, 1e+2]);
xticks([1e-3 1e-2 1e-1 1 10 100]);
ylim([0.335 1]);
yticks([0 0.4 0.5 0.6 0.7 0.8 0.9 1]);
yticklabels({' ','0.4','0.5','0.6','0.7','0.8','0.9','1'});

xlabel('r (distance)');
set(gca, 'XScale', 'log');
title('Numerical solution VS Analytical Approximation');
subtitle('\lambda = 0.1,\nu = 10, \eta = 1');
legend('Numerical', 'Padé2','Exp-Var','EXP-Global','DblExp-Var','DblExp-Global','Location', 'southeast');
