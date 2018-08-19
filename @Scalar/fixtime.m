function evalObj = fixtime(obj, tau)
%FIXTIME - partial evaluation of Scalar along the first (time) dimension.
%
%   evalObj = FIXTIME(obj, tau) returns a Scalar of dimension (d-1) corresponding to evaluation of obj along its first dimension. In usual cases
%       this dimension is utilized by the time variable so that this evaluation is the time-tau map for the Scalar.
%       IMPORTANT NOTE: This evaluation is in MATERIAL time, not real time.
%
%   Inputs:
%       obj - Scalar
%       tau - double
%
%   Outputs:
%       evalObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 23-Apr-2017; Last revision: 13-Aug-2018

%   TODO:
%   support for Fourier and Chebyshev bases
%   support for evaluation along other dimensions
%   support for specifying real time evaluation

if ~strcmp(obj.Basis, 'Taylor')
    warning('fixtime - only implemented for Taylor basis')
end


if length(obj) ==1 % obj is a single Scalar
    switch obj.Dimension
        case 2
            time_vector = bsxfun(@power,tau,0:size(obj.Coefficient,1)-1);
            evalObj = Scalar(time_vector*obj.Coefficient, obj.Basis);
        case 3
            if tau==1
                evalObj = Scalar(squeeze(sum(obj.Coefficient,1)), obj.Basis); % sum over the first dimension is equivalent to evaluation at tau = 1.
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

end
end % fixtime

% Revision History:
%{
15-May-2017 - support for 2-dimensional Scalars
13-Aug-2018 - updated for Scalar class
%}
