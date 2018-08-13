function convObj = fouriertaylortimes(obj, rightObj)
%FOURIERTAYLORTIMES - mixed convolution for Fourier-Taylor series
%
%   Syntax:
%       output = FOURIERTAYLORTIMES(input1, input2)
%       output = FOURIERTAYLORTIMES(input1, input2, input3)
%
%   Description:
%       FOURIERTAYLORTIMES() - description
%
%   Inputs:
%       input1 - Description
%       input2 - Description
%       input3 - Description
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Aug-2018; Last revision: 04-Aug-2018

%%

switch obj.Dimension
    case 2
        [M,N] = size(obj.Coefficient);
        if ~isequal([M,N], size(rightObj.Coefficient))
            error('mixed basis convolution not tested for non-uniform sized Scalars')
        end

        if strcmp(obj.NumericalClass, 'double')
            linearConv = conv2(obj.Coefficient, rightObj.Coefficient); % full linear convolution in both directions
            circularConv = linearConv(1:M,1:N); % truncate in the Taylor (time) direction
            circularConv(:,1:N-1) = circularConv(:,1:N-1) + linearConv(1:M, N+1:end); % circular convolution in Fourier (spatial) direction
            convObj = Scalar(circularConv);
        elseif strcmp(obj.NumericalClass, 'intval')

        end
    otherwise
        error('not implemented yet')

end
end % end fouriertaylortimes

