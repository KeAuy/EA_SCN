%APPROXIMATIONS Evaluate a set of closed-form approximations of b(r).
%
%   APPROX = BUFFER.APPROXIMATIONS(R, LAMBDA, NU, ETA, NAMES, OPTIONS)
%   returns one field per entry of NAMES; every field is NUMEL(R)-by-NUMEL(NU).
%   NAMES may be Pade, Pade2, RBA2, Linear, ExpSer, ExpVar, ExpGlobal,
%   ExpPade, DblExpSer, DblExpVar and DblExpGlobal.
%
%   OPTIONS: PadeRefine, PadeRefineMethod ('fsolve' or 'newton'),
%   PadeRealCoefficients and RealBranchIndex.

function approx = approximations(r, lambda, nu, eta, names, options)

arguments
    r double {mustBePositive}
    lambda (1,1) double {mustBePositive}
    nu double
    eta (1,1) double
    names = {}
    options.PadeRefine (1,1) logical = false
    options.PadeRefineMethod (1,:) char ...
        {mustBeMember(options.PadeRefineMethod, {'fsolve', 'newton'})} = 'fsolve'
    options.PadeRealCoefficients (1,1) logical = true
    options.RealBranchIndex (1,1) double {mustBeNonnegative} = 0
end

if isempty(names)
    names = {'Pade', 'Pade2', 'RBA2', 'Linear', 'ExpSer', 'ExpVar', ...
             'ExpGlobal', 'ExpPade', 'DblExpSer', 'DblExpVar', ...
             'DblExpGlobal'};
end
if ischar(names) || isstring(names)
    names = cellstr(names);
end

r  = r(:);
nu = nu(:);
q  = 1 ./ (eta + nu.');

approx = struct();

for k = 1:numel(names)
    name = names{k};
    switch lower(name)
        case 'pade'
            approx.pade = buffer.pade1(r, lambda, nu, eta);
        case 'pade2'
            approx.pade2 = buffer.pade2(r, lambda, nu, eta, ...
                'Refine', options.PadeRefine, ...
                'RefineMethod', options.PadeRefineMethod, ...
                'RealCoefficients', options.PadeRealCoefficients);
        case 'rba2'
            approx.rba2 = buffer.rba2(r, lambda, nu, eta);
        case 'linear'
            approx.linear = buffer.expSeries(r, lambda, nu, eta, 0);
        case 'expser'
            approx.expSeries = buffer.expSeries(r, lambda, nu, eta, 1 / 2);
        case 'expvar'
            S = (1 + 2 * q * eta) / 3;
            approx.expVar = buffer.expSeries(r, lambda, nu, eta, S);
        case 'expglobal'
            S = log(3 / 2) + q * eta * log(4 / 3);
            approx.expGlobal = buffer.expSeries(r, lambda, nu, eta, S);
        case 'exppade'
            approx.expPade = buffer.expPade(r, lambda, nu, eta);
        case 'dblexpser'
            approx.doubleExpSeries = buffer.doubleExpSeries( ...
                r, lambda, nu, eta, options.RealBranchIndex);
        case 'dblexpvar'
            approx.doubleExpVar = buffer.doubleExpVariant( ...
                r, lambda, nu, eta, 'var');
        case 'dblexpglobal'
            approx.doubleExpGlobal = buffer.doubleExpVariant( ...
                r, lambda, nu, eta, 'global');
        otherwise
            error('buffer:unknownApproximation', ...
                'Unknown approximation name: %s.', name);
    end
end
end
