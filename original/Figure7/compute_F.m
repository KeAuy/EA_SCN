function F = compute_F_Fig6(b0, b, r, hL, hR, lambda, nu, eta)
%COMPUTE_F_Fig6  compute_F_Fig2 的向量化版本（对 lambda 同时算）。
%   b      : N x M   每一列是一个 lambda 的解
%   lambda : 1 x M   （标量也行，会自动广播）
%   b0     : 1 x M 或标量
%   F      : N x M
%
%   与标量版逐元素相同，只是把 k 的循环换成了矩阵运算。

[N, M] = size(b);
if isscalar(lambda), lambda = lambda*ones(1, M); end
if isscalar(b0),     b0     = b0*ones(1, M);     end

q = 1/(eta + nu);

F = zeros(N, M);
F(1,:) = b(1,:) - b0;

k   = (2:N-1).';
hb  = 0.5*(hL(k) + hR(k));

F(k,:) = lambda .* ((b(k+1,:) - b(k,:))./hR(k) - (b(k,:) - b(k-1,:))./hL(k))./hb ...
       + 2*lambda .* (b(k+1,:) - b(k-1,:))./((hL(k) + hR(k)).*r(k)) ...
       - (b(k,:) - 1).*(nu*b(k,:) + eta) - b(k,:)./r(k);

F(N,:) = b(N,:) - (1 - q/r(N) + eta*q^3/r(N)^2);
end