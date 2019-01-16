function coefDecay = decay(obj)
%DECAY - returns the vector of coefficient norms in the first dimension
%
%   coefDecay = DECAY(oby) returns the ell_1 norm of each coefficient of the Scalar with respect to the first dimension.
%
%   Inputs:
%       obj - Scalar
%
%   Outputs:
%       coefDecay - double
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 09-Jun-2017; Last revision: 13-Aug-2018

%   TODO:
%   Add support for Fourier and Chebyshev coefficients.

if numel(obj) > 1 % vectorized method
    rowDecay = arrayfun(@(j)obj(j).decay, 1:numel(obj));
    coefDecay = reshape(rowDecay, size(obj));
    
elseif strcmp(obj.NumericalClass, 'Scalar')
    coefDecay = obj.Coefficient(end).norm();
    
else % single Scalar
    coefDecay = sumdim(mid(obj.Coefficient), 2:obj.Dimension);
end
end % decay

function sumDim = sumdim(coef, dimension)
% sum across each dimension specified
if isequal(length(dimension), 1)
    sumDim = sum(abs(coef), dimension);
else
    sumDim = sumdim(sum(abs(coef),dimension(end)), dimension(1:end-1)); % sum last dimension and recurse
end
end

% Revision History:
%{
02-Jul-2017 - support for Scalar coefficient type
13-Aug-2018 - updated for Scalar class
21-Aug-2018 - support for arbitrary dimensions
%}

