function objTailRatio = tailratio(obj,tailBeginsAt)
% computes the norm of the higher order terms of a polynomial relative to the polynomial's norm

% Written by S. Kepley 04/2017
% Added support for intvals 08/2017

% INPUT:
% tailBeginsAt (double): vector of values specifying where the tail terms begin

objDimsidx = obj(1).dims();
Sz = size(obj(1).Coefficient);
Sz = Sz(objDimsidx);
if any(tailBeginsAt < 1) % Non-tail terms specified as fraction of total terms.
    fractionTerms = tailBeginsAt < 1;
    tailBeginsAt(fractionTerms) = round(tailBeginsAt(fractionTerms).*Sz(fractionTerms));
    tailBeginsAt(tailBeginsAt < 1) = ones(size(tailBeginsAt < 1));
end

if length(obj) > 1 % vectors of Scalars
%     getRatio = @(j)obj(j).tailratio(tailBeginsAt);
    objTailRatio(length(obj)) = obj(end).tailratio(tailBeginsAt);
    for j = 1:length(obj)-1
       objTailRatio(j) = obj(j).tailratio(tailBeginsAt);
    end
else
    switch length(tailBeginsAt)
        case 1 % obj is 1d
            objTail = Scalar(obj.Coefficient(tailBeginsAt:end));

        case 2 % obj is 2d
            tailCoefficient = obj.Coefficient;
            tailCoefficient(1:tailBeginsAt(1)-1, 1:tailBeginsAt(2)-1) = zeros(1:tailBeginsAt(1)-1, 1:tailBeginsAt(2)-1);
            objTail = Scalar(tailCoefficient);

        case 3 % obj is 3d
            objTail = Scalar(obj.Coefficient(tailBeginsAt(1):end,tailBeginsAt(2):end,tailBeginsAt(3):end));

        otherwise
            error('objTailRatio not implemented for Scalar with dimension greater than 3')
    end
    tailNorm = objTail.norm();
    objTailRatio = tailNorm./obj.norm;
end
if isa(objTailRatio,'intval')
    objTailRatio = sup(objTailRatio);
end
end
