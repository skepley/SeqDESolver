function convObj = intvaltimes(leftObj, rightObj)
%INTVALTIMES - Scalar multiplication (linear convolution) for intval coefficients

%   Syntax:
%       output = INTVALTIMES(input1, input2)
%       output = INTVALTIMES(input1, input2, input3)
%
%   Description:
%       INTVALTIMES() - description
%
%   Inputs:
%       leftObj - Scalar objects with intval coefficients
%       rightObj - Scalar objects with intval coefficients
%
%   Outputs:
%       convObj - Scalar corresponding to convolution leftObj*rightObj
%
%   Subfunctions: none
%   Classes required: @Scalar
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 25-Jul-2018; Last revision: 25-Jul-2018

%%
switch leftObj.Dimension
    case 2
        [leftCoefficient, leftExponent] = leftObj.exponent();
        [rightCoefficient, rightExponent] = rightObj.exponent();

        leftExpArray = repmat(leftExponent,[1 1 length(rightCoefficient)]); %1-by-1-by-#Coefficients in right factor
        rightExpArray = shiftdim(repmat(rightExponent',[1 1 length(leftCoefficient)]),2); %1-by-1-by-#Coefficients in left factor

        productExponent = 1 + shiftdim(leftExpArray + rightExpArray,2);
        coefArray = leftCoefficient*rightCoefficient'; % outer product
        convObj = full(sparse(productExponent(:,:,1), productExponent(:,:,2), coefArray));
    otherwise
        error('not yet implemented')

end
end % end intvaltimes

