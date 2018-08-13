function coefDecay = decay(obj)
% returns norm of last Taylor coefficient collapsed onto first variable.


% Written S. Kepley 06/2017
% Updated for Scalar coefficient type 07/2017

% TODO: Add support for measuring decay in other dimensions.

% ---------------------- INPUT ----------------------
% obj is a Scalar
% varargin is dimension in which to measure the decay

% ---------------------- OUTPUT ----------------------
% decay of the coefficients in the required dimension

if numel(obj) > 1 % obj is a vector of Scalars
    rowDecay = arrayfun(@(j)obj(j).decay,1:numel(obj));
    coefDecay = reshape(rowDecay,size(obj));
elseif strcmp(obj.NumericalClass,'Scalar')
    coefDecay = obj.Coefficient(end).norm();
else
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
                finalCoefficient = Scalar(obj.Coefficient(end,:),obj.Truncation(2:end));
                coefDecay = finalCoefficient.norm();
            else
                finalCoefficient = Scalar(obj.Coefficient(end,:),obj.Truncation(2:end),obj.Weight(2:end));
                coefDecay = finalCoefficient.norm();
            end
        case 3
            if strcmp(obj.Weight,'ones')
                finalCoefficient = Scalar(obj.Coefficient(end,:,:),obj.Truncation(2:end));
                coefDecay = finalCoefficient.norm();
            else
                finalCoefficient = Scalar(obj.Coefficient(end,:,:),obj.Truncation(2:end),obj.Weight(2:end));
                coefDecay = finalCoefficient.norm();
            end
    end
    if isa(obj.Coefficient,'intval')
        coefDecay = mid(coefDecay);
    end
end
end
