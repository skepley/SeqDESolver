function scalarNorm = norm(obj, varargin)
%NORM - computes the weighted ell_1 norm of a Scalar
%
%   NORM returns the weighted ell_1 norm for a Scalar along specified dimensions. The weights are specified in the Scalar properties.
%
%   Syntax:
%       scalarNorm = NORM(obj) returns the ell_1 norm of obj
%       scalarNorm = NORM(obj, dims) returns the ell_1 norm of obj in the directions specified by dims
%
%   Inputs:
%       obj - Scalar
%       dims - Vector of dimensions along which to sum
%
%   Outputs:
%       scalarNorm - double or intval depending on Scalar NumericalClass
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

if numel(obj) > 1 % vectorized norm
    for j = 1:numel(obj)
        scalarNorm(j) = obj(j).norm(varargin{:});
    end
    scalarNorm = reshape(scalarNorm,size(obj));

else
    if strcmp(obj.Weight,'ones')
        if nargin == 2
            normDim = varargin{1};
            scalarNorm = sum(abs(obj.Coefficient),normDim);
        else
            scalarNorm = sum(abs(obj.Coefficient(:)));
        end

    elseif obj.Dimension ~= 2
        error('norm not implemented for this weight and dimension')

    else
        weight_matrix = bsxfun(@(x,y)obj.Weight(1).^(x).*obj.Weight(2).^(y),0:obj.Truncation(2)-1,(0:obj.Truncation(1)-1)');
        scalarNorm = sum(dot(weight_matrix,abs(obj.Coefficient)));
    end
end
end % end norm

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
%}
