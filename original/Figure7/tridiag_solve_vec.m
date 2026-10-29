function x = tridiag_solve_vec(aL, aD, aR, f)
%TRIDIAG_SOLVE_VEC  一次解 M 个三对角系统（Thomas 算法，全列并行）。
%   aL,aD,aR : N x M  下/主/上对角线
%   f        : N x M  右端
%   x        : N x M  解
%   aL(1,:) 和 aR(N,:) 不参与运算。

[N, M] = size(aD);

cp = zeros(N, M);
dp = zeros(N, M);
cp(1,:) = aR(1,:)./aD(1,:);
dp(1,:) = f(1,:)./aD(1,:);

for k = 2:N
    m = aD(k,:) - aL(k,:).*cp(k-1,:);
    cp(k,:) = aR(k,:)./m;
    dp(k,:) = (f(k,:) - aL(k,:).*dp(k-1,:))./m;
end

x = zeros(N, M);
x(N,:) = dp(N,:);
for k = N-1:-1:1
    x(k,:) = dp(k,:) - cp(k,:).*x(k+1,:);
end
end