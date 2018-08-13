function scaletime(obj,L)
% rescales time

% Written by S. Kepley 05/2017
% Updated for Scalar coefs 07/2017

% ---------------------- INPUT ----------------------
% obj (Scalar): A Taylor series with 1st coordinate (time) scaled to t = 1
% L (double): The new time rescaling

% ---------------------- OUTPUT ----------------------
% obj (Scalar): The rescaled object

if numel(obj) > 1 % vectorized norm
    %                 scale_matr = Scalar(repmat(bsxfun(@power,abs(L),[0:obj(1).Truncation(1)-1]'),[1,obj(1).Truncation(2)]));
    for j = 1:length(obj)
        obj(j).scaletime(L);
    end
elseif obj.Truncation(1) ==1
    error('obj.Truncation should not equal 1') % rule out constant with respect to time
elseif isa(obj.Coefficient,'Scalar')
    for k = 1:obj.Truncation(1)
        obj.Coefficient(k).scaletime(L^(k-1));
    end
else
    switch obj.Dimension
        case 1
            obj.Coefficient = abs(L)*obj.Coefficient;
        case 2
            obj.Coefficient = repmat(bsxfun(@power,abs(L),[0:obj.Truncation(1)-1]'),[1,obj.Truncation(2)]).*obj.Coefficient;
        case 3
            scaleCoefficient = @(j)L^(j-1).*obj.Coefficient(j,:,:);
            for j = 1:obj.Truncation(1)
                obj.Coefficient(j,:,:) = scaleCoefficient(j);
            end
        otherwise
            error('Not implemented for higher dimension')
    end
end
end
