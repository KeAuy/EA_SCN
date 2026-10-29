function res = bcfun(ya, yb, lambda)
%BCFUN  Boundary residuals of the buffer-gas BVP, for use with BVP5C.
%
%   RES = BCFUN(YA, YB, LAMBDA) evaluates the two boundary conditions of the
%   radial problem described in ODEFUN:
%
%       b'(rMin) = b(rMin)/(2*lambda)   % near-origin behaviour
%       b(rMax)  = 1                    % undisturbed buffer far away
%
%   YA and YB are the state vectors [b; db/dr] at the inner and outer end
%   points, respectively.
%
%   See also ODEFUN, GUESS, SHOOTING_FIG2.

res = [ya(2) - ya(1)/(2*lambda);
       yb(1) - 1];

end
