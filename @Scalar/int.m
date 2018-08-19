function intObj = int(obj, varargin)
%INT - shift Scalar coefficients by one in the first dimension.
%
%   Syntax:
%       intObj = INT(obj) returns a Scalar which has been "integrated" in the first dimension. On the coefficient level, the m^th coefficient
%       has been shifted and rescaled by 1/m. This is embedded into the truncation space which has 1 additional coefficient and the first
%       coefficient is padded by zero.
%
%   Inputs:
%       obj - Scalar
%
%   Outputs:
%       intObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

if ~strcmp(obj.Basis, 'Taylor')
    error('int - non-Taylor bases not implemented')
end

if nargin > 1
    nIntegration = varargin{1}; % iterate operator n times
else
    nIntegration = 1;
end

if nIntegration > 1
    obj = int(obj, nIntegration-1);
end

switch obj.Dimension
    case 1        
        switch obj.NumericalClass
            case 'intval'
                scaleCoefficient = objCoefficient./intval(1:obj.Truncation(1))';
                shiftCoefficient = [intval(0); scaleCoefficient];
            case 'double'
                scaleCoefficient = obj.Coefficient./(1:obj.Truncation(1))';
                shiftCoefficient = [0; scaleCoefficient];
            otherwise
                error('int')
        end
        
    case 2
        intObj = Scalar([zeros(1,obj.Truncation(2));obj.Coefficient(1:obj.Truncation(1)-1,:)], obj.Basis);
        error('int - this has not been updated for 2-d yet')
        
    otherwise
        error('shift method not implemented in higher dimensions')
end
intObj = Scalar(shiftCoefficient, obj.Basis);
end %  shift

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
18-Aug-2018 - changed name from "shift" to "int", added coefficient rescaling to the old shifting operation
%}
