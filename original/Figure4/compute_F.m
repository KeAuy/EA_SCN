% function F = compute_F_Fig2(b0,b, r,h_minus,h_plus, N,lambda,nu,eta)
% q = 1/(eta+nu);
% F = zeros(N,1);
% F(1) = b(1)-b0;
% for i = 2:N-1
%     ri = r(i);
%     hm = h_minus(i);
%     hp = h_plus(i);
% 
%     brr = (2/(hm+hp))*...
%           ((b(i+1)-b(i))/hp ...
%           -(b(i)-b(i-1))/hm);
% 
%     br = ( hm^2*(b(i+1)-b(i)) ...
%          + hp^2*(b(i)-b(i-1)) ) ...
%          /(hp*hm*(hp+hm));
% 
%     F(i) = -lambda*(brr+(2/ri)*br) ...
%            +(b(i)-1)*(nu*b(i)+eta) ...
%            +b(i)/ri;
% end
% 
% 
% F(N) = b(N) - (1 - (q/r(end)) + (eta*(q.^3))/ (r(end)^2) ...
%     + (eta*(1-2*q*eta).*(q.^4))/(r(end)^3) + (eta*(2*lambda + 5*eta*(q.^2)*(eta*q- 1)*(q.^4)))/(r(end)^4));
% end


function F = compute_F_Fig2(b0, b, r, hL, hR, N, lambda, nu, eta)
% Residual of  lambda*(b'' + 2 b'/r) = (b-1)(nu*b+eta) + b/r
% on the non-uniform (geometric) grid r(1:N), with hL(k)=r(k)-r(k-1),
% hR(k)=r(k+1)-r(k).  Second-order accurate in s = ln r.
%
% Boundary conditions:
%   b(1) = b0                        (value at the inner cut-off)
%   b(N) = 1 - q/r(N) + eta*q^3/r(N)^2   (two-term far-field asymptotics, Eq.22)
%
% NOTE the far-field condition: imposing b(N)=1 at a finite rmax is only
% first-order accurate and lifts the whole error curve (it kills the dip of
% Fig. 2C).  If you prefer b(N)=1, use rmax >= 1e6.

q = 1/(eta+nu);
F = zeros(N,1);
F(1) = b(1) - b0;

for k = 2:N-1
    hb = 0.5*(hL(k)+hR(k));
    F(k) = lambda*((b(k+1)-b(k))/hR(k) - (b(k)-b(k-1))/hL(k))/hb ...
         + 2*lambda*(b(k+1)-b(k-1))/((hL(k)+hR(k))*r(k)) ...
         - (b(k)-1)*(nu*b(k)+eta) - b(k)/r(k);
end

F(N) = b(N) - (1 - q/r(N) + eta*q^3/r(N)^2);
end