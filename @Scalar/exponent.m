function [coefficient, exponent] = exponent(obj)
%EXPONENT - Return a vector of coefficients and exponents for specified Scalar
%
%   Syntax:
%       output = EXPONENT(input1, input2)
%       output = EXPONENT(input1, input2, input3)
%
%   Description:
%       EXPONENT() - description
%
%   Inputs:
%       obj - instance of Scalar object
%
%   Outputs:
%       coefficient - coefficient array as a column vector
%       exponent - exponent array of size: #length(coefficients)-by-Scalar dimension
%
%   Subfunctions: none
%   Classes required: @Scalar
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Aug-2018; Last revision: 04-Aug-2018

%%
coefficient = reshape(obj.Coefficient,[],1); % coefficient as column vector

% exponent arrays
switch obj.Dimension
    case 1
        exponent = (0:obj.Truncation(1))';

    case 2
        [var1Exp,var2Exp] = meshgrid(0:obj.Truncation(1) -1,0:obj.Truncation(2) -1);
        exponent = [reshape(var2Exp,[],1), reshape(var1Exp,[],1)];

    otherwise
        error('not implemented yet')
end
end % end exponent

