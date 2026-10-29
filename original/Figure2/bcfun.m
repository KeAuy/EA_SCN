function res = bcfun(ya, yb, lambda)

res = [ya(2) - ya(1)/(2*lambda);
       yb(1) - 1];

end
