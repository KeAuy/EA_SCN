
function res = bcfun(ya, yb, lambda, nu)

eta = 10;                
q   = 1/(eta + nu);
rR  = 1e2;               

res = [ ya(2) - ya(1)/(2*lambda);
        yb(1) - (1 - q/rR + eta*q^3/rR^2) ];
end
