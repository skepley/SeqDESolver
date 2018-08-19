function varargout = sqrt(obj, varargin)
%SQRT - compute square roots for Scalars using automatic differentiation.
%
%   B = SQRT(obj) is a scalar with the same truncation embedding as obj which approximately satisfies B*B = obj. The constant term for B is positive.
%
%   B = SQRT(obj, -1) is the negative solution branch (i.e. B*B = obj and the constant term for B is negative).
%
%   B = SQRT(obj) is a Scalar with the same truncation embedding as obj which approximately satisfies B*B = obj.
%
%   [B, inv(B)] is a pair of Scalars satisfying B*B = obj and B*inv(B) = 1.
%
%   Inputs:
%       obj - Scalar
%       branch - integer either {1, -1}
%
%   Outputs:
%       sqrtObj - Scalar
%       invObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-May-2018; Last revision: 13-Aug-2018

if ~strcmp(obj.Basis, 'Taylor')
    warning('sqrt - only implemented for Taylor basis')
end


if strcmp(obj.NumericalClass,'intval')
    warning('sqrt - untested for interval coefficients')
end


if abs(obj.Coefficient(1,1)) < 1e-9
    warning('constant term < 1e-9. square root should not be trusted')
end


if nargin ==2
    branch = varargin{1};
else
    branch = 1;
end

switch obj.Dimension
    case 0 % sqrt of a real scalar
        a0 = abs(obj.Coefficient(1));
        sqrtObj = Scalar(branch*sqrt(a0), obj.Basis);
        invObj = Scalar(1/sqrtObj.Coefficient(1), obj.Basis);

    case 1 % sqrt of 1d power series
        dObj = obj.ds;
        a0 = abs(obj.Coefficient(1));
        sqrtObj = Scalar(branch*sqrt(a0), obj.Basis, 1,1);
        invObj = Scalar(1/sqrtObj.Coefficient(1),1,1);
        for m = 1:obj.Truncation - 1
            dA = Scalar(dObj.Coefficient(1:m), obj.Basis, invObj.Truncation);
            dsqrt = dA*invObj; % obj'*
            uu = invObj*invObj;
            sqrtNext = .5*(1/m)*dsqrt.Coefficient(m);
            invNext = .5*(-1/m)*mtimes(uu,dsqrt,'Recursion');
            sqrtObj.append(sqrtNext);
            invObj.append(invNext);
        end
    case 2 % sqrt of 2d power series
        dObj = obj.dt;
        dObj.Coefficient = dObj.Coefficient(2:end,:);
        a0 = Scalar(obj.Coefficient(1,:) obj.Basis);
        [sqrtInit,invInit] = sqrt(a0);
        sqrtObj = Scalar(sqrtInit, obj.Basis, [1,a0.Truncation]);
        invObj = Scalar(invInit, obj.Basis, [1,a0.Truncation]);
        for m = 1:obj.Truncation(1)-1
            dA = Scalar(dObj.Coefficient(1:m,:), obj.Basis, invObj.Truncation);
            dsqrt = dA*invObj;
            uu = invObj*invObj;
            sqrtNext = .5*(1/m)*dsqrt.Coefficient(m,:);
            invNext = .5*(-1/m)*mtimes(uu,dsqrt,'Recursion');
            sqrtObj.append(sqrtNext);
            invObj.append(invNext);
        end
    otherwise
        error('sqrt not implemented for this dimension')
end

if nargout == 1
    varargout{1} = sqrtObj;
elseif nargout ==2
    varargout{1} = sqrtObj;
    varargout{2} = invObj;
end
end

% Revision History:
%{
17-Jun-2017 - computation rewritten to use automatic differentiation
    support for lower dimensional square roots
    added option to return inv(obj) which is computed as a by product of automatic differentiation
13-Aug-2018 - updated for Scalar class
%}
