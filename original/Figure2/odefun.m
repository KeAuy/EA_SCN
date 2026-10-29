function dbdr = odefun(r, y, nu, lambda, eta)
%ODEFUN  Right-hand side of the radial buffer-gas equation.
%
%   The second-order problem
%
%       lambda*(b'' + 2*b'/r) = (b - 1)*(nu*b + eta) + b/r
%
%   is integrated as the first-order system y = [b; db/dr]:
%
%       db/dr   = db
%       ddb/dr  = -(2/r)*db + ((b - 1)*(nu*b + eta) + b/r)/lambda
%
%   See also BCFUN, GUESS.

b  = y(1);
db = y(2);

dbdr = [db;
        -(2/r)*db + ((b - 1)*(nu*b + eta) + b/r)/lambda];

end
