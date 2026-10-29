%MEANERROR Radial mean of the pointwise error between two profiles.
%
%   MODE 'abs' gives mean(|approx - reference|); 'log' gives
%   mean(|log(approx) - log(reference)|).  The result is a column vector.

function err = meanError(reference, approx, mode)

arguments
    reference double
    approx double
    mode (1,:) char {mustBeMember(mode, {'abs', 'log'})} = 'abs'
end

if strcmpi(mode, 'log')
    pointwise = abs(log(approx) - log(reference));
else
    pointwise = abs(approx - reference);
end

err = mean(pointwise, 1).';
end
