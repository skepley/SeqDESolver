function invObj = inv(obj,varargin)
% Compute b = 1/A for specified A by solving the differential equation u' = -A'u^2

% TODO: Add validation and make the argument recursive.

% Written by S. Kepley 05/2017
% Updated surfaceDimension 06/2017
% Support for Scalar coef 07/2017

if strcmp(obj.NumericalClass,'Scalar')
    const = obj.Coefficient(1).Coefficient(1);
else
    const = obj.Coefficient(1);
end

if abs(const) < 1e-13 % check if constant term is near 0.
    disp('Scalar is not invertible. Inverse should not be trusted')
end

switch obj.Dimension
    case 0 % inverse of a real scalar
        a0 = obj.Coefficient(1);
        invObj = Scalar(1/a0);

    case 1 % inverse of 1d power series
        ddtObj = obj.ds;
        a0 = obj.Coefficient(1);
        invObj = Scalar(1/a0,obj.Truncation);
        for m = 1:obj.Truncation - 1
            uu = invObj*invObj;
            Aprime = ddtObj.Coefficient(m:-1:1);
            newCoefficient = -(1/m)*dot(uu.Coefficient(1:m),Aprime);
            invObj.Coefficient(m+1) = newCoefficient;
        end

    case 2 % inverse of 2d power series
        deg = size(obj.Coefficient);
        ddtObj = obj.dt;
        ddtObj.Coefficient = [ddtObj.Coefficient(2:end,:);zeros(1,deg(2))]; % This padding shouldn't be necessary.
        a0 = Scalar(obj.Coefficient(1,:));
        invObj = Scalar(inv(a0),[1,a0.Truncation]);
        for m = 1:obj.Truncation(1)-1
            uu = invObj*invObj;
            Aprime = Scalar(ddtObj.Coefficient(1:m,:),uu.Truncation);
            newCoefficient = -(1/m)*mtimes(uu,Aprime,'Recursion');
            invObj.append(newCoefficient);
        end
    otherwise
        error('not yet implemented')
end
end

