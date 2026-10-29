function [aL, aD, aR] = compute_J_Fig6(b, r, hL, hR, lambda, nu, eta)
%COMPUTE_JDIAG_FIG2_VEC  compute_J_Fig2 的三条对角线（向量化）。
%   不再组装稀疏矩阵，直接给 [下对角 aL, 主对角 aD, 上对角 aR]，每条 N x M。
%   之后用 tridiag_solve_vec 一次解掉 M 个三对角系统。
%
%   aL(1,:) 与 aR(N,:) 不参与运算，置 0 即可。

[N, M] = size(b);
if isscalar(lambda), lambda = lambda*ones(1, M); end

aL = zeros(N, M);
aD = ones(N, M);      % 第一行、最后一行是 Dirichlet
aR = zeros(N, M);

k  = (2:N-1).';
hb = 0.5*(hL(k) + hR(k));

aL(k,:) = lambda./hL(k)./hb - 2*lambda./((hL(k) + hR(k)).*r(k));
aR(k,:) = lambda./hR(k)./hb + 2*lambda./((hL(k) + hR(k)).*r(k));
aD(k,:) = (-lambda./hR(k) - lambda./hL(k))./hb ...
        - (2*nu*b(k,:) - nu + eta) - 1./r(k);
end