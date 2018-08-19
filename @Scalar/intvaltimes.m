function convObj = intvaltimes(obj, rightObj)
%INTVALTIMES - Scalar multiplication (linear convolution) for intval coefficients
%
%   Description:
%       convObj = INTVALTIMES(obj, rightObj) - description
%
%   Inputs:
%       obj - Scalar object with intval coefficients
%       rightObj - Scalar object with intval coefficients
%
%   Outputs:
%       convObj - Scalar corresponding to linear convolution obj*rightObj
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 25-Jul-2018; Last revision: 18-Aug-2018

%   TODO:
%   support for higher dimensions
%   support for other bases

if ~strcmp(obj.Basis, rightObj.Basis)
    error('intvaltimes - Scalar multiplication requires same basis')
end

if ~strcmp(obj.Basis, 'Taylor')
    error('intvaltimes - this multiplication is only implemented for Taylor coefficients')
end

switch obj.Dimension
    case 1 % 1-d discrete convolution
        [leftCoefficient, leftExponent] = obj.exponent();
        [rightCoefficient, rightExponent] = rightObj.exponent();
        leftExpArray = repmat(leftExponent, [1 length(rightCoefficient)]); 
        rightExpArray = repmat(rightExponent',[length(leftCoefficient),1]); 
        
        productExponent = 1 + leftExpArray+rightExpArray;
        coefArray = leftCoefficient*rightCoefficient'; % outer product
        convObjCoefficient = full(sparse(productExponent, 1, coefArray));
        
    case 2 % 2-d discrete convolution
        [leftCoefficient, leftExponent] = obj.exponent();
        [rightCoefficient, rightExponent] = rightObj.exponent();
        leftExpArray = repmat(leftExponent, [1, 1, length(rightCoefficient)]); %1-by-1-by-#Coefficients in right factor
        rightExpArray = shiftdim(repmat(rightExponent', [1, 1, length(leftCoefficient)]),2); %1-by-1-by-#Coefficients in left factor
        
        productExponent = 1 + shiftdim(leftExpArray + rightExpArray,2);
        coefArray = leftCoefficient*rightCoefficient'; % outer product
        convObjCoefficient = full(sparse(productExponent(:,:,1), productExponent(:,:,2), coefArray));
    otherwise
        error('not yet implemented')
        
end
convObj = Scalar(convObjCoefficient, obj.Basis);

end % intvaltimes

% Revision History:
%{
13-Aug-2018 - updated for Scalar class
18-Aug-2018 - support for 1-d convolution

%}
