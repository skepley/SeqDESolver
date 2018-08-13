function dObj_dt = dt(obj)
%DT - compute the derivative of a Scalar
%
%   DT(obj) returns a Scalar corresponding to the derivative of obj in the first dimension. The derivative is embedded in the same
%   truncation space by padding with zeros.
%
%   Inputs:
%       obj - Scalar object
%
%   Outputs:
%       dObj_dt - Scalar object
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

if obj.Dimension ==1 % space is 0-dim
    C = 0:obj.Truncation(1)-1;
    dObj_dt = Scalar(C.*obj.Coefficient);

elseif obj.Dimension == 2 % space is 1-dim
    C = repmat((0:obj.Truncation(1)-1),obj.Truncation(2),1)';
    dObj_dt = Scalar(C.*obj.Coefficient);

elseif obj.Dimension ==3 % space is 2-dim
    C = zeros(obj.Truncation);
    for j = 2:obj.Truncation(1)
        C(j,:,:) = (j-1)*ones(obj.Truncation(2:3));
    end
    dObj_dt = Scalar(C.*obj.Coefficient);
end
end % end dt

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
%}
