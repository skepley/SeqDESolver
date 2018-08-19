function coefDecay = decay(obj)
% returns norm of last coefficient collapsed onto first variable.
%
%   coefDecay = DECAY(oby) returns the ell_1 norm of the last coefficient in Scalar with respect to the first variable.
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
%   Add support for measuring decay in other dimensions.
%   Add support for Fourier and Chebyshev coefficients.

if numel(obj) > 1 % vectorized method
    rowDecay = arrayfun(@(j)obj(j).decay, 1:numel(obj));
    coefDecay = reshape(rowDecay,size(obj));

elseif strcmp(obj.NumericalClass, 'Scalar')
    coefDecay = obj.Coefficient(end).norm();

else % single Scalar
    switch obj.Dimension
        case 0
            coefDecay = abs(obj.Coefficient(1));
        case 1
            if strcmp(obj.Weight,'ones')
                coefDecay = abs(obj.Coefficient(end));
            else
                coefDecay = obj.Weight*abs(obj.Coefficient(end));
            end
        case 2
            if strcmp(obj.Weight,'ones')
                finalCoefficient = Scalar(obj.Coefficient(end,:), obj.Basis, obj.Truncation(2:end));
                coefDecay = finalCoefficient.norm();
            else
                finalCoefficient = Scalar(obj.Coefficient(end,:), obj.Basis, obj.Truncation(2:end),obj.Weight(2:end));
                coefDecay = finalCoefficient.norm();
            end
        case 3
            if strcmp(obj.Weight,'ones')
                finalCoefficient = Scalar(obj.Coefficient(end,:,:), obj.Basis, obj.Truncation(2:end));
                coefDecay = finalCoefficient.norm();
            else
                finalCoefficient = Scalar(obj.Coefficient(end,:,:), obj.Basis, obj.Truncation(2:end), obj.Weight(2:end));
                coefDecay = finalCoefficient.norm();
            end
    end % switch

    if isa(obj.Coefficient,'intval')
        coefDecay = mid(coefDecay);
    end
end
end % decay

% Revision History:
%{
02-Jul-2017 - support for Scalar coefficient type
13-Aug-2018 - updated for Scalar class
%}

