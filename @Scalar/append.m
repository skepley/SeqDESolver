function append(obj,nextCoefficient,varargin)
% Appends a new coefficient to given dimension (default is 1st dimension).

% Written by S. Kepley 05/9/17
% Updated for multiple dimensions 06/24/17
% Updated for Scalar coefficients 07/02/17

% ---------------------- INPUT ----------------------
% obj is a Scalar
% nextCoefficient is a coefficient array of size [1,obj.Dimension-1]

% ---------------------- OUTPUT ----------------------
% obj (Scalar): The rescaled object

if strcmp(obj.NumericalClass,'Scalar') % append Scalar coefficient
    obj.Coefficient(end+1) = nextCoefficient;
    obj.Truncation(1) = obj.Truncation(1) + 1;

else % append double/intval array
    if isa(nextCoefficient,'Scalar')
        nextCoefficient = nextCoefficient.Coefficient;
    end

    try
        switch obj.Dimension
            case 1
                obj.Coefficient(end+1) = nextCoefficient;
            case 2
                obj.Coefficient(end+1,:) = nextCoefficient;
            case 3
                obj.Coefficient(end+1,:,:) = nextCoefficient;
            otherwise
                error('Not yet implemented for higher dimension')
        end
        obj.Truncation(1) = obj.Truncation(1) + 1;
    catch
        error('Appended coefficient must have same dimension as the surface being appended to')
    end
end
end
