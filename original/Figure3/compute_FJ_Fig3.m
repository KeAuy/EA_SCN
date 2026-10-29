%% 辅助函数：计算残差和雅可比（修正边界条件）
function [F, J] = compute_FJ(b, r, h, N, lambda, nu, eta, b0)
    q = 1/(eta+nu);
    F = zeros(N,1);
    J = sparse(N,N);
    
    % r=0 边界条件
    F(1) = b(1) - b0;
    J(1,1) = 1;
    
    % 内点
    for i = 2:N-1
        ri = r(i);
        brr = (b(i+1) - 2*b(i) + b(i-1)) / h^2;
        br  = (b(i+1) - b(i-1)) / (2*h);
        F(i) = -lambda*(brr + (2/ri)*br) + (b(i)-1)*(nu*b(i)+eta) + b(i)/ri;
        
        % 雅可比（基于中心差分）
        J(i, i-1) = -lambda * (1/h^2 - 1/(ri*h));
        J(i, i)   =  (2*lambda)/h^2 + 2*nu*b(i) + eta - nu + 1/ri;
        J(i, i+1) = -lambda * (1/h^2 + 1/(ri*h));
    end
    
    % 无穷远边界条件
    %F(N) = b(N) - 1;
    F(N) = b(N) - (1 - (q/r(end)) + (eta*(q.^3))/ (r(end)^2) ...
    + (eta*(1-2*q*eta).*(q.^4))/(r(end)^3) + (eta*(2*lambda + 5*eta*(q.^2)*(eta*q- 1)*(q.^4)))/(r(end)^4));
    J(N,N) = 1;
end