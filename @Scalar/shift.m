function shiftObj = shift(obj)
%SHIFT - shift Scalar coefficients by one in the first dimension.
%
%   Syntax:
%       shiftObj = SHIFT(obj) returns a Scalar whose coefficients have been shifted. For Taylor series this corresponds to multiplication by
%       the first variable. This is embedded into the same truncation space so the last coefficient is dropped and the first is padded by zero.
%
%   Inputs:
%       obj - Scalar
%
%   Outputs:
%       shiftObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

switch obj.Dimension
    case 1
        shiftObj = Scalar([0, obj.Coefficient(1:end-1)]);

    case 2
        shiftObj = Scalar([zeros(1,obj.Truncation(2));obj.Coefficient(1:obj.Truncation(1)-1,:)]);

    otherwise
        error('shift method not implemented in higher dimensions')
end
end % end shift

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
%}
