%INITIALPROFILE Piecewise initial profile for the Newton iterations.
%
%   [B, SWITCHRADIUS] = BUFFER.INITIALPROFILE(R, B0, LAMBDA, NU, ETA, SWITCH)
%   uses the small-r expansion below SWITCH and the far-field expansion
%   elsewhere.  An empty SWITCH picks the radius where the two expansions are
%   relatively closest.

function [b, switchRadius] = initialProfile(r, b0, lambda, nu, eta, switchRadius)

arguments
    r double {mustBePositive}
    b0 double
    lambda (1,1) double {mustBePositive}
    nu double
    eta (1,1) double
    switchRadius double = []
end

r  = r(:);
nu = nu(:)';
[bNear, ~, bFar, ~] = buffer.expansion(r, b0, lambda, nu, eta, 4);

if isempty(switchRadius)
    gap = abs(bNear - bFar) ./ (abs(bNear) + abs(bFar) + 1e-30);
    [~, idx] = min(gap, [], 1);
    switchRadius = r(idx).';
else
    if isscalar(switchRadius)
        switchRadius = repmat(switchRadius, 1, numel(nu));
    else
        switchRadius = switchRadius(:).';
    end
end

mask = r < switchRadius;
b    = bFar;
b(mask) = bNear(mask);
end
