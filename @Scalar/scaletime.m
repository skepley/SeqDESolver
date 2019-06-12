function scaledObj = scaletime(obj, L)
%SCALETIME - rescale domain of Scalar in the first (time) dimension by z --> Lz.
%
%   scaledObj = SCALETIME(obj, L) returns a Scalar which defines a function f_L: D^d ---> R which is equivalent to the
%       function defined by obj, f:[-L,L] x D^{d-1} ---> R.
%
%   Inputs:
%       obj - Scalar
%       L - field scalar (double or intval)
%
%   Outputs:
%       scaledObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-May-2017; Last revision: 13-Aug-2018

if ~strcmp(obj(1).Basis, 'Taylor')
    warning('scaletime - only implemented for Taylor basis')
end


% if strcmp(obj(1).NumericalClass, 'intval')
%     warning('scaletime - rescaling interval coefficients may cause excessive wrapping.')
% end

if numel(obj) > 1 % vectorized norm
    scaledObj = Scalar;
    for j = 1:length(obj)
        scaledObj(j) = obj(j).scaletime(L);
    end
    scaledObj = reshape(scaledObj, size(obj));
elseif obj.Truncation(1) ==1
    error('obj.Truncation should not equal 1') % rule out constant with respect to time
else
    switch obj.Dimension
        case 1
            scaledObj = Scalar(abs(L)*obj.Coefficient, obj.Basis);
        case 2
            scaledObj = Scalar(repmat(bsxfun(@power,abs(L),(0:obj.Truncation(1)-1)'),[1,obj.Truncation(2)]).*obj.Coefficient, obj.Basis);
        case 3
            error('Not fixed for dimensions 3 and above')
            scaleCoefficient = @(j)L^(j-1).*obj.Coefficient(j,:,:);
            for j = 1:obj.Truncation(1)
                obj.Coefficient(j,:,:) = scaleCoefficient(j);
            end
        otherwise
            error('Not implemented for higher dimension')
    end
end
end

% Revision History:
%{
02-Jul-2018 - support for Scalar coefficients
13-Aug-2018 - updated for Scalar class
%}
