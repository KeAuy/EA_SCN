clear; clc;

eta     = 1;
b0_test = 0.5;
N       = 100;
Nlam    = 100;
Nnu     = 100;

r        = 10.^(-3 + 5*(1:N)/N);
lambda_v = logspace(-3, 2, Nlam);
nu_v     = logspace(-3, 2, Nnu);

options = bvpset('NMax', 50000, 'RelTol', 1e-6, 'AbsTol', 1e-6);

b0 = zeros(Nlam, Nnu);

use_pb = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
if use_pb
    dq      = parallel.pool.DataQueue;
    t_start = tic;
    afterEach(dq, makeProgress(Nlam, lambda_v, t_start));
end

parfor i = 1:Nlam
    lam     = lambda_v(i);
    row     = zeros(1, Nnu);
    solprev = [];

    for k = 1:Nnu
        nu = nu_v(k);

        if isempty(solprev)
            solinit = bvpinit(r, @(rr) guess(b0_test, lam, nu, eta, rr));
        else
            solinit = solprev;
        end

        solprev = bvp5c(@(rr, y) odefun(rr, y, nu, lam, eta), ...
                        @(ya, yb) bcfun(ya, yb, lam), ...
                        solinit, options);

        Y      = deval(solprev, r);
        row(k) = Y(1, 1);
    end

    b0(i, :) = row;
    if use_pb, send(dq, i); end
end

save('b0_Fig6.mat', 'b0');
