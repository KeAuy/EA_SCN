function err = meanError(reference, approx, mode)
%MEANERROR Radial mean of the pointwise error between two profiles.
%
%   ERR = BUFFER.MEANERROR(REFERENCE, APPROX, MODE) returns the column
%   vector of r-averaged errors.  MODE is 'abs' (default) for
%   mean(|approx-reference|) or 'log' for
%   mean(|log(approx)-log(reference)|).
%
%   See also BUFFER.PLOTERRORS.

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
