% function J = compute_J_Fig2(b,r,h_minus,h_plus,N,lambda,nu,eta)
% 
% J = sparse(N,N);
% 
% J(1,1) = 1;
% J(N,N) = 1;
% 
% for i = 2:N-1
% 
%     hm = h_minus(i);
%     hp = h_plus(i);
%     ri = r(i);
% 
%     % ===== b'' coefficients =====
%     A =  2/(hm*(hm+hp));   % b_{i-1}
%     B = -2/(hm*hp);        % WRONG form ignored (not used)
%     C =  2/(hp*(hm+hp));   % b_{i+1}
% 
%     % more stable: recompute correct central coefficient
%     B = -(A + C);  % ensures consistency
% 
%     % ===== b' coefficients (EXACT from your F) =====
%     D = -hp^2/(hp*hm*(hp+hm));   % b_{i-1}
%     E = (hp^2 - hm^2)/(hp*hm*(hp+hm)); % b_i
%     G =  hm^2/(hp*hm*(hp+hm));   % b_{i+1}
% 
%     % ===== linear operator part =====
%     J(i,i-1) = -lambda*(A + (2/ri)*D);
% 
%     J(i,i) = -lambda*(B + (2/ri)*E) ...
%              + (2*nu*b(i) + eta) ...
%              + 1/ri ...
%              - nu;   % from derivative of (b-1)(nu b + eta)
% 
%     J(i,i+1) = -lambda*(C + (2/ri)*G);
% 
% end
% end


function J = compute_J_Fig2(b, r, hL, hR, N, lambda, nu, eta)

J = spalloc(N, N, 3*N);
J(1,1) = 1;
J(N,N) = 1;

for k = 2:N-1
    hb = 0.5*(hL(k)+hR(k));
    J(k,k-1) = lambda/hL(k)/hb - 2*lambda/((hL(k)+hR(k))*r(k));
    J(k,k)   = (-lambda/hR(k) - lambda/hL(k))/hb ...
               - (2*nu*b(k) - nu + eta) - 1/r(k);
    J(k,k+1) = lambda/hR(k)/hb + 2*lambda/((hL(k)+hR(k))*r(k));
end
end
