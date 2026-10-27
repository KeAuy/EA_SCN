function y0 = buffer_guess(r, b0, lambda, nu, eta)
%BUFFER_GUESS Initial guess for BVP5C.
%
%   Y0 = BUFFER_GUESS(R, B0, LAMBDA, NU, ETA) returns the state [b; db/dr] at
%   the radii in R. The small-r expansion is used where the two expansions of
%   BUFFER_EXPANSION agree to within 10%, and the large-r expansion
%   elsewhere. Near the switch-over radius the small-r branch is preferred.
%
%   R may be a scalar or a vector; Y0 is 2-by-NUMEL(R).
%
%   See also BUFFER_EXPANSION, BUFFER_ODE, BUFFER_BC.

[bInner, dbInner, bOuter, dbOuter] = buffer_expansion(r, b0, lambda, nu, eta);

agreement = abs(bInner - bOuter) ./ (abs(bInner) + abs(bOuter) + 1e-14);
useInner  = agreement < 0.1 | r < 1;

b  = bOuter;
db = dbOuter;
b(useInner)  = bInner(useInner);
db(useInner) = dbInner(useInner);

y0 = [b(:).'; db(:).'];
end
