function dbdr = odefun(r, y, nu, lambda, eta)

b  = y(1);
Db = y(2);

dbdr = zeros(2,1);

dbdr(1) = Db;

dbdr(2) = -(2/r)*Db ...
          + (1/lambda)*(b-1)*(nu*b + eta) ...
          + b/(lambda*r);
end