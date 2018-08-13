function evalObj = fixtime(obj,tau)
% collapse Taylor series onto fixed time t

% Written S. Kepley 04/2017
% Added 2D funcionality 05/15/17

% ---------------------- INPUT ----------------------


% ---------------------- OUTPUT ----------------------

if length(obj) ==1 % obj is a single Scalar
    switch obj.Dimension
        case 2
            time_vector = bsxfun(@power,tau,0:size(obj.Coefficient,1)-1);
            evalObj = Scalar(time_vector*obj.Coefficient);
        case 3
            if tau==1
                evalObj = Scalar(squeeze(sum(obj.Coefficient,1))); % sum over the first dimension is equivalent to evaluation at tau = 1.
            elseif tau == -1

            else
                error('Not implemented for arbitrary tau in 2D')
            end
    end

else % obj is a vector of Scalars
    evalObj(length(obj)) = obj(length(obj)).fixtime(tau);
    for j = 1:length(obj)-1
        evalObj(j) = obj(j).fixtime(tau);
    end

%     newCoefficient = arrayfun(@(j)obj(j).fixtime(tau),1:length(obj),'UniformOutput',false);
%     evalObj = newCoefficient{1};
%     for j = 2:length(newCoefficient)
%         if obj(j).Dimension == 1
%             evalObj(j) = newCoefficient{j};
%         else
%             evalObj = cat(1,evalObj,newCoefficient{j});
%         end
%     end
%
%     if obj(1).Dimension == 1
%         evalObj = reshape(evalObj,size(obj)); % output shape matches input shape
%     end
end
end
