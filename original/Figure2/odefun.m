function dbdr = odefun(r, y, nu, lambda, eta)

b  = y(1);
db = y(2);

dbdr = [db;
        -(2/r)*db + ((b - 1)*(nu*b + eta) + b/r)/lambda];

end
