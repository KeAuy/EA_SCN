clear,clc;
eta = 1;
N  = 800;
N1 = 100;
lambda_vec = logspace(-3,2,N1)';
nu_vec = logspace(-3, 2, N1)';

tol      = 1e-3;
max_iter = 500;

k_idx  = (1:N)';
r_vec  = 10.^(-3 + 5.*k_idx./N);
q_vec  = 1 ./ (eta + nu_vec);
x_vec  = 1 ./ r_vec;

alpha_r = 10^(5/N);
h_L = r_vec .* (1 - 1/alpha_r);
h_R = r_vec .* (alpha_r - 1);

NLAM = numel(lambda_vec);
NNU  = numel(nu_vec);
LAM  = reshape(lambda_vec, 1, NLAM);

b_num = zeros(N,NLAM, NNU);
res_l = zeros(NLAM, NNU);

temp  = load("b0_Fig6.mat");
b_r0  = temp.b0;

for j = 1:NNU
    nu  = nu_vec(j);
    q   = 1/(eta + nu);
    b0j = reshape(b_r0(:,j), 1, NLAM);

    c1v = b0j./(2*LAM);
    c2v = ((b0j - 1).*(nu*b0j + eta) + c1v)./(6*LAM);
    b_frob = b0j + c1v.*r_vec + c2v.*r_vec.^2;

    xv = 1./r_vec;
    b_asymp = 1 - q*xv + eta*q^3*xv.^2 + eta*(1 - 2*eta*q)*q^4*xv.^3 ...
            + eta*(2*LAM + 5*eta*q^2*(eta*q - 1))*q^4.*xv.^4;

    rel = abs(b_frob - b_asymp)./(abs(b_frob) + abs(b_asymp));
    [~, isw] = min(rel, [], 1);
    rsw = reshape(r_vec(isw), 1, NLAM);
    msk = r_vec < rsw;
    b   = double(msk).*b_frob + double(~msk).*b_asymp;

    xtol = 1e-12;

    for it = 1:max_iter
        F   = compute_F(b0j, b, r_vec, h_L, h_R, LAM, nu, eta);
        res = max(abs(F), [], 1);
        if max(res) < tol, break; end

        [aL, aD, aR] = compute_J(b, r_vec, h_L, h_R, LAM, nu, eta);
        delta = tridiag_solve_vec(aL, aD, aR, -F);

        stp = max(abs(delta), [], 1) ./ (1 + max(abs(b), [], 1));
        if all(stp < xtol), break; end

        t = ones(1, NLAM);
        for ls = 1:max_iter
            bt  = b + t.*delta;
            Ft  = compute_F(b0j, bt, r_vec, h_L, h_R, LAM, nu, eta);
            rt  = max(abs(Ft), [], 1);
            bad = (stp >= xtol) & (rt > (1 - 1e-4*t).*res);
            if ~any(bad), break; end
            t(bad) = 0.5*t(bad);
        end
        b = b + t.*delta;
    end

    b_num(:,:,j) = b;
    F = compute_F(b0j, b, r_vec, h_L, h_R, LAM, nu, eta);
    res_l(:,j) = max(abs(F), [], 1);
    fprintf('nu = %.3e 完成：Newton %2d 步，max|F| = %.2e，最大相对步长 %.1e\n', ...
            nu, it, max(res_l(:,j)), max(stp));
end
save('b_num_Fig6_A1A2.mat',"b_num");

LAM = reshape(lambda_vec(:), 1, [], 1);
NU  = reshape(nu_vec(:),     1, 1, []);
QV  = 1 ./ (eta + NU);
PV  = QV .* NU;

Pc_R = 24*LAM - 3.*PV.*QV.*(PV - 2);
Pc_K = 8*LAM.^2 ...
     - 2*LAM.*QV.*(2*PV.^2 - 5*PV - 2) ...
     - PV.*QV.^2.*(PV - 2);
Pc_Q = Pc_K + 24*LAM.^2 - 4*LAM.*PV.*QV.*(2*PV - 3);
Pc_V = 2*LAM.*(PV + 6) + PV.*QV.*(PV - 3).*(PV - 2);
Pc_F = 6*LAM.^2 ...
     - LAM.*QV.*(PV + 3).*(PV - 2) ...
     - eta.*PV.*QV.^3.*(2*PV - 3);
Pc_J = 24*LAM.^2 ...
     - 2*LAM.*PV.*QV.*(5*PV - 8) ...
     - eta.*(PV.^2).*QV.^3.*(PV - 2);
Pc_E = 12*LAM.^2 + eta.*PV.*QV.^2.*(8*LAM + eta.*PV.*QV.^2);
Pc_H = (PV.^2).*(QV.^4).*((PV - 2).^2) + 1600*LAM.^4 ...
     - 104*(LAM.^3).*QV.*(PV - 2).*(9*PV + 4) ...
     + 4*(LAM.^2).*(QV.^2).*(39*PV.^4 - 87*PV.^3 - 63*PV.^2 + 184*PV + 4) ...
     - 4*LAM.*PV.*(QV.^3).*(PV - 2).*(6*PV.^3 - 21*PV.^2 + 20*PV + 2);
Pc_W = 12000*LAM.^5 ...
     - 4*(LAM.^4).*QV.*(2687*PV.^2 - 2966*PV - 553) ...
     + 4*(LAM.^3).*(QV.^2).*(729*PV.^4 - 1626*PV.^3 + 486*PV.^2 + 421*PV + 70) ...
     - (LAM.^2).*(QV.^3).*(243*PV.^6 - 810*PV.^5 + 609*PV.^4 + 210*PV.^3 ...
        - 21*PV.^2 - 240*PV - 4) ...
     - 2*LAM*eta.*PV.*(QV.^5).*(PV - 2).*(20*PV.^2 - 13*PV - 1) ...
     + (PV.^2).*(eta^2).*(QV.^7);

Pc_G = Pc_Q.^3 + 2*LAM.*Pc_R.*(3*Pc_Q.*Pc_F - 9*LAM.*Pc_E.*Pc_R ...
     + Pc_V.*sqrt(complex(-3.*QV.*Pc_W)));
cu   = Pc_G.^(1/3);

B2 = real((QV./Pc_R).*(Pc_Q + cu + Pc_H./cu));
B1 = (1./Pc_V).*(((B2.*(B2.*Pc_R - 3.*QV.*Pc_K))./(6*LAM.*QV.^2)) - Pc_J);

for itP = 1:40
    fP1 = 2*LAM.*B1.^2 + B1.*B2 - B2.^2./QV ...
        - 2*LAM.*QV.^2.*B1 - (QV.^2 + 2*LAM).*B2;
    fP2 = (1 + QV).*B1.^2 - 2*B1.*B2./QV + 2*(3*LAM - QV.^3).*B1 ...
        + (1 - QV).*B2 - QV.^2.*(6*LAM + QV.^2.*(1 - QV));
    resP = max(max(abs(fP1(:))), max(abs(fP2(:))));
    if ~isfinite(resP) || resP < 1e-10, break; end

    jP11 = 4*LAM.*B1 + B2 - 2*LAM.*QV.^2;
    jP12 = B1 - 2*B2./QV - (QV.^2 + 2*LAM);
    jP21 = 2*(1 + QV).*B1 - 2*B2./QV + 2*(3*LAM - QV.^3);
    jP22 = -2*B1./QV + (1 - QV);
    detP = jP11.*jP22 - jP12.*jP21;
    B1 = B1 + (-fP1.*jP22 + fP2.*jP12)./detP;
    B2 = B2 + (-fP2.*jP11 + fP1.*jP21)./detP;
end

A1 = B1 - QV;
A2 = B2 - QV.*(B1 - eta*QV.^2);

RV = reshape(r_vec(:), [], 1, 1);
b_Pade2 = (RV.^2 + A1.*RV + A2) ./ (RV.^2 + B1.*RV + B2);

RV  = reshape(r_vec(:),       [], 1, 1);
LAM = reshape(lambda_vec(:),  1, [], 1);
NU  = reshape(nu_vec(:),      1, 1, []);
QV  = 1 ./ (eta + NU);

D = (1 + QV./RV).^2 - 4*NU.*QV.^2./RV;

b_RBA = 1 - (1./(2*NU.*QV)) .* (1 + QV./RV - sqrt(D));
b_RBA = b_RBA + 0*LAM;

b_RBA2 = b_RBA + 2*LAM*eta ./ ((1 + RV./QV).^2 - 4*NU.*RV).^2;

RV  = reshape(r_vec(:),      [], 1, 1);
LAM = reshape(lambda_vec(:), 1, [], 1);
NU  = reshape(nu_vec(:),     1, 1, []);
QV  = 1 ./ (eta + NU);

P_exp_pade = 6 + (eta*(QV.^2).*(1 - eta*QV))./LAM;
Q_exp_pade = (1./P_exp_pade).*(6 + (eta*(QV.^2))./(2*LAM));
R_exp_pade = Q_exp_pade.*(1 - 2*LAM./QV) + (2*eta^2*(QV.^2))./P_exp_pade;
U_exp_pade = R_exp_pade/3 - (Q_exp_pade.^2)/2;
V_exp_pade = Q_exp_pade.*(R_exp_pade - Q_exp_pade.^2) + (12*LAM)./(QV.*P_exp_pade);
W_exp_pade = (Q_exp_pade.^2).*(R_exp_pade/3 - Q_exp_pade.^2/4) ...
           + 4*LAM.*((6*(Q_exp_pade + LAM./QV) - eta^2*QV.^2)./(3*QV.*P_exp_pade));

E_exp_pade = (V_exp_pade.^2).*((V_exp_pade.^2)/4 ...
           + U_exp_pade.*(U_exp_pade.^2 - 3*W_exp_pade)) ...
           - W_exp_pade.*(W_exp_pade - 3*U_exp_pade.^2).^2;

G_exp_pade = V_exp_pade.^2/2 + U_exp_pade.*(U_exp_pade.^2 - 3*W_exp_pade) ...
           + sqrt(complex(E_exp_pade));
cG = G_exp_pade.^(1/3);

H_exp_pade = (W_exp_pade + U_exp_pade.^2)./cG + cG - 2*U_exp_pade;

alpha_exp_pade = (sqrt(H_exp_pade) - Q_exp_pade ...
    + sqrt(2*V_exp_pade./sqrt(H_exp_pade) - H_exp_pade - 6*U_exp_pade))./(2*LAM);

beta_exp_pade = eta*QV.^3 ./ (QV.*alpha_exp_pade.*(LAM.*alpha_exp_pade + 1) - 1);

b_exp_pade = 1 + QV.*((exp(-alpha_exp_pade.*RV) - 1)./RV) ...
               + (QV.^3*eta)./(beta_exp_pade + RV.^2);

RV  = reshape(r_vec(:),      [], 1, 1);
LAM = reshape(lambda_vec(:), 1, [], 1);
NU  = reshape(nu_vec(:),     1, 1, []);
QV  = 1 ./ (eta + NU);

P_Dbl_Exp_Var = LAM*(8*log(2) - 5) ...
              + ((4/3)*(QV.^2)*eta).*(1 - eta*QV).*(1 - 3*log(4/3));
Q_Dbl_Exp_Var = LAM + ((2/3)*(QV.^2)*eta) ...
              .* (1 - 6*log(9/8) + (2*QV*eta).*(1 - 6*log(4/3)));
R_Dbl_Exp_Var = (QV*eta + 2)/3;

W_Dbl_Exp_Var = 3*eta*(QV.^2).*P_Dbl_Exp_Var;
E_Dbl_Exp_Var = (Q_Dbl_Exp_Var.^2) + W_Dbl_Exp_Var.*R_Dbl_Exp_Var;
F_Dbl_Exp_Var = (3/2)*(Q_Dbl_Exp_Var.*R_Dbl_Exp_Var ...
              - W_Dbl_Exp_Var./QV).*W_Dbl_Exp_Var + Q_Dbl_Exp_Var.^3;

pk   = [1/2 + (sqrt(3)/2)*1i;  1/2 - (sqrt(3)/2)*1i;  -1];
Cbrt = (sqrt(complex(F_Dbl_Exp_Var.^2 - E_Dbl_Exp_Var.^3)) ...
        + F_Dbl_Exp_Var).^(1/3);

AL3 = zeros(3, numel(lambda_vec), numel(nu_vec));
for k = 1:3
    Gk = pk(k) * Cbrt;
    AL3(k,:,:) = reshape((Q_Dbl_Exp_Var - Gk - E_Dbl_Exp_Var./Gk) ...
                         ./ W_Dbl_Exp_Var, 1, numel(lambda_vec), numel(nu_vec));
end

ALr   = reshape(AL3, 3, []);
b0_3d = 1 - QV.*AL3 + 0.5*eta*QV.^3.*AL3.^2;
b0r   = reshape(b0_3d, 3, []);

ok = abs(imag(ALr)) < 1e-10*max(1,abs(ALr)) & real(ALr) > 0 ...
     & real(b0r) > 0 & real(b0r) < 1;

score = real(ALr);  score(~ok) = Inf;
[best1, kb] = min(score, [], 1);

score2 = real(ALr); score2(real(ALr) <= 0) = Inf;
[best2, kb2] = min(score2, [], 1);
need1 = ~isfinite(best1);  kb(need1) = kb2(need1);

need2 = ~isfinite(best1) & ~isfinite(best2);
if any(need2), [~, kb3] = max(real(ALr), [], 1); kb(need2) = kb3(need2); end

ALPHA_Dbl_Exp_Var = reshape(ALr(sub2ind(size(ALr), kb, 1:size(ALr,2))), ...
                            1, numel(lambda_vec), numel(nu_vec));

DblExpVar = 1 + QV.*((exp(-ALPHA_Dbl_Exp_Var.*RV) - 1)./RV) ...
          - (QV.^3)*eta.*((exp(-ALPHA_Dbl_Exp_Var.*RV) ...
             .*(1 + ALPHA_Dbl_Exp_Var.*RV) - 1)./RV.^2);

RV  = reshape(r_vec(:),      [], 1, 1);
LAM = reshape(lambda_vec(:), 1, [], 1);
NU  = reshape(nu_vec(:),     1, 1, []);
QV  = 1 ./ (eta + NU);

P_Dbl_Exp_Global = 2*LAM*(1 - log(2)) + QV.^2*eta.*(1 - QV*eta)*(log(3) - 1);
Q_Dbl_Exp_Global = LAM - 2*QV.^2*eta.*(1 - log(81/32) + 2*QV*eta*log(9/8));
R_Dbl_Exp_Global = QV*eta + 2*(1 - QV*eta)*log(3/2);

W_Dbl_Exp_Global = 3*eta*QV.^2.*P_Dbl_Exp_Global;
E_Dbl_Exp_Global = Q_Dbl_Exp_Global.^2 + W_Dbl_Exp_Global.*R_Dbl_Exp_Global;
F_Dbl_Exp_Global = (3/2)*(Q_Dbl_Exp_Global.*R_Dbl_Exp_Global ...
                 - W_Dbl_Exp_Global./QV).*W_Dbl_Exp_Global ...
                 + Q_Dbl_Exp_Global.^3;

pk1 = 1/2 + (sqrt(3)/2)*1i;
G_Dbl_Exp_Global = pk1 * (sqrt(complex(F_Dbl_Exp_Global.^2 ...
                    - E_Dbl_Exp_Global.^3)) + F_Dbl_Exp_Global).^(1/3);
ALPHA_Dbl_Exp_Global = (Q_Dbl_Exp_Global - G_Dbl_Exp_Global ...
                     - E_Dbl_Exp_Global./G_Dbl_Exp_Global) ./ W_Dbl_Exp_Global;

DblExpGlobal = 1 + QV.*((exp(-ALPHA_Dbl_Exp_Global.*RV) - 1)./RV) ...
             - (QV.^3)*eta.*((exp(-ALPHA_Dbl_Exp_Global.*RV).*(1 + ALPHA_Dbl_Exp_Global.*RV) - 1)./RV.^2);

RV  = reshape(r_vec(:),      [], 1, 1);
LAM = reshape(lambda_vec(:), 1, [], 1);
NU  = reshape(nu_vec(:),     1, 1, []);
QV  = 1 ./ (eta + NU);

S_expglobal     = log(3/2) + QV*eta*log(4/3);
alpha_expglobal = (sqrt(S_expglobal.^2 + LAM./QV) - S_expglobal)./LAM;

ExpGlobal = 1 + QV.*((exp(-alpha_expglobal.*RV) - 1)./RV);

methods = {'Pade2', 'RBA2','Exp-Pade',  'Dbl-ExpGlobal','ExpGlobal'};
B       = {b_Pade2, b_RBA2, b_exp_pade,  DblExpGlobal, ExpGlobal};

nM   = numel(methods);
Nlam = numel(lambda_vec);
Nnu  = numel(nu_vec);

E = nan(Nlam, Nnu, nM);

for i = 1:nM
    A = B{i};
    D = abs(real(A) - b_num);
    E(:,:,i) = reshape(mean(D, 1), Nlam, Nnu);
end

[Emin, winner] = min(E, [], 3);

methods = {'Pade2', 'RBA2', 'Exp-Pade','DblExp-Global','Exp-Global'};
nameList = {'Pade2', 'RBA2','Exp-Pade',  'DblExp-Global','Exp-Global'};

colList  = [ 0.13 0.69 0.30;
             0.55 0.10 0.10;
             0.50 0.50 0.50;
             0.9961 0.4844  0.9922;
             0.00 0.00 0.90];

nM   = numel(methods);
cmap = zeros(nM, 3);
for i = 1:nM
    k = find(strcmp(nameList, methods{i}), 1);
    if isempty(k)
        warning('方法 "%s" 没有指定颜色，用黑色代替', methods{i});
        cmap(i,:) = [0 0 0];
    else
        cmap(i,:) = colList(k,:);
    end
end

NaNcol = [0.97 0.97 0.97];
cmapFull = [cmap; NaNcol];

W = winner;

figure('Color','w','Position',[100 100 720 560]);
imagesc(nu_vec, lambda_vec, W);
set(gca,'XScale','log','YScale','log');
axis xy;
clim([0.5, nM + 1.5]);
colormap(gca, cmapFull);
cb = colorbar; cb.Ticks = 1:nM+1; cb.TickLabels = [methods, {'none'}];
xlabel('\nu'); ylabel('\lambda');
title(sprintf('最优方法分区  (\\eta = %g)', eta));
colormap(gca, cmapFull);
cb = colorbar;
cb.Ticks      = 1:nM+1;
cb.TickLabels = [methods, {'none'}];
cb.FontSize   = 9;
box on;
exportgraphics(gcf, 'winner_map.png', 'Resolution', 160);
fprintf('已保存 winner_map.png\n');

A_v = eta * (QV.^2) .* P_Dbl_Exp_Var;
B_v = -Q_Dbl_Exp_Var;
C_v = -R_Dbl_Exp_Var;
D_v = 1 ./ QV;

Delta_v = 18*A_v.*B_v.*C_v.*D_v - 4*B_v.^3.*D_v + B_v.^2.*C_v.^2 - 4*A_v.*C_v.^3 - 27*A_v.^2.*D_v.^2;

Delta_v = reshape(Delta_v, Nlam, Nnu);
s_v = max(abs(Delta_v(:)));
if s_v > 0
    Delta_v = Delta_v / s_v;
end

A_g = eta * (QV.^2) .* P_Dbl_Exp_Global;
B_g = -Q_Dbl_Exp_Global;
C_g = -R_Dbl_Exp_Global;
D_g = 1 ./ QV;

Delta_g = 18*A_g.*B_g.*C_g.*D_g - 4*B_g.^3.*D_g + B_g.^2.*C_g.^2 - 4*A_g.*C_g.^3 - 27*A_g.^2.*D_g.^2;

Delta_g = reshape(Delta_g, Nlam, Nnu);
s_g = max(abs(Delta_g(:)));
if s_g > 0
    Delta_g = Delta_g / s_g;
end

ax = gca;
hold(ax, 'on');

if any(Delta_v(:) < 0) && any(Delta_v(:) > 0)
    contour(ax, nu_vec, lambda_vec, Delta_v, [0 0], ...
            'LineColor', [0.85 0.85 0.85], 'LineWidth', 1.0);
end

if any(Delta_g(:) < 0) && any(Delta_g(:) > 0)
    contour(ax, nu_vec, lambda_vec, Delta_g, [0 0], ...
            'LineColor', [0.85 0.65 0.00], 'LineWidth', 1.0);
end
hold(ax, 'off');

Z = log10(Emin);

figure('Color','w');
contourf(nu_vec, lambda_vec, Z);

set(gca, 'XScale','log', 'YScale','log', ...
         'XLim',[1e-3 1e2], 'YLim',[1e-3 1e2], ...
         'XTick',10.^(-3:2), 'YTick',10.^(-3:2), ...
         'XTickLabel',{'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'}, ...
         'YTickLabel',{'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});

colormap(bone); clim('auto');

cb = colorbar; cb.Ticks = -6:1:-2;

xlabel('\nu'); ylabel('\lambda');
title('log_{10}||b_{num}-b_{approx}||, \eta=1');
