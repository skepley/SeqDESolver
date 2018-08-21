function productObj = mtimes(leftObj,rightObj,varargin)
%MTIMES - Define multiplication of Scalars
%
%   Syntax:
%   C = MTIMES(obj, B) is the Scalar convolution product when B is a Scalar truncated to size of obj.
%
%   C = MTIMES(T, obj) is the Scalar image of T(obj) when T is an Operator (linear operator over sequence space) truncated to size of obj.
%
%   C = MTIMES(a, obj) is the Scalar given by a*obj when a is a double truncated to size of obj.
%
%   C = MTIMES(obj, a) is the Scalar given by a*obj = obj*a when a is a double or intval (intvals require right multiplication) truncated to size of obj.
%
%   C = MTIMES(A,B, 'Full') returns a Scalar with no truncation.
%
%   C = MTIMES(obj, B, 'Recursion') returns on the (m-1)^st coefficient with respect to the first dimension when obj.Truncation(1) = m.
%
%   C = MTIMES(A,B, truncationArray) returns a Scalar truncated in each dimension by truncation Array = [N1,...,Nd]).
%
%   Inputs:
%       leftObj - Scalar, double, or Operator
%       rightObj - Scalar, double, or intval
%       truncationSize - string (Full or Recursion) or integer vector of length d.
%
%   Outputs:
%       productObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Mar-2017; Last revision: 13-Aug-2018


if isa(leftObj,'Operator') % Operator acts on Scalar (left only)
    if ~strcmp(rightObj.Basis, 'Taylor')
        warning('mtimes - only implemented for Taylor basis')
    end

    switch size(obj,1)
        case 1
            Lv = leftObj.Matrix*rightObj.col;
            productObj = reshape(Lv,rightObj.Truncation);
        otherwise
            if size(leftObj) ~= size(obj,1)*ones(1,2)
                error('Operator matrix and vector of Scalars must have the same size')
            else
                % productCell = arrayfun(@(j)dot(leftObj(j,:),obj),1:size(obj,1),'UniformOutput',false);
                % productObj = cell2mat(productCell);
            end
    end

elseif isa(leftObj,'Scalar') && isa(rightObj,'Scalar') % Multiplication of 2 Scalars (Cauchy product of analytic functions)
    if ~strcmp(rightObj.Basis, 'Taylor')
        warning('mtimes - only implemented for Taylor basis')
    end

    if ~strcmp(leftObj.Basis, 'Taylor')
        warning('mtimes - only implemented for Taylor basis')
    end
    
    if ~strcmp(leftObj.Basis, rightObj.Basis)
        error('mtimes - Scalars must have the same basis type')
    else
       basis = leftObj.Basis; 
    end

    % this is currently a catch all that needs to be fixed
    if isequal(leftObj.Truncation,[]) || isequal(rightObj.Truncation,[])
        disp('Scalar should never have empty value for Truncation')
        L = Scalar(leftObj.Coefficient, basis);
        R = Scalar(rightObj.Coefficient, basis);
        productObj = Scalar(mtimes(L, R, varargin{:}), basis);
        productObj.Truncation = [];
        return
    % one factor is a constant function
    elseif (leftObj.Dimension == 0 || rightObj.Dimension == 0)
        productObj = Scalar(leftObj.Coefficient*rightObj.Coefficient, basis);
        return
    % products require factors of the same dimension
    elseif ~isequal(leftObj.Dimension, rightObj.Dimension)
        error('Cauchy products for factors of different dimensions is not defined yet')

    % left/right Dimension match. Set product Dimension equal
    else
        Dimension = leftObj.Dimension;
    end

    if (nargin > 2) % If both factors have double NumericalClass default behavior is to truncate the product to the same size as obj. Otherwise call with varargin.
        truncationSize = varargin{1}; % truncation types: {'Fixed','Recursion','Full',ArraySize}
        % 'Fixed' requires two Scalars with double NumericalClass of identical size.
        % 'Recursion' gives fast convolution for computing Taylor coefficient by recursion. Returns only the (M+1)-st coefficient in 1st dimension. Note: Object returned in this case is a double, not a Scalar.
        % 'Full' produces the full Cauchy product including higher order.
        % 'FFT' same as Fixed but uses FFT for the product.
        % ArraySize = [M,N1,N2,...] produces products truncated to the size specified.

    elseif ~isa(leftObj.Coefficient,'intval') && ~isa(rightObj.Coefficient,'intval') && ~isequal(leftObj.Truncation,rightObj.Truncation)
        disp('Scalars arent the same size. Using Full option.')
        truncationSize = 'Full';
    else
        truncationSize = 'Fixed';
    end

    switch Dimension
        case 1 % ---------------------------------------- 1D SURFACES ----------------------------------
            if strcmp(leftObj.NumericalClass,'intval') || strcmp(rightObj.NumericalClass,'intval')
                if strcmp(truncationSize,'Recursion')
                    minMode = min([leftObj.Truncation,rightObj.Truncation]);
                    productObj = dot(leftObj.Coefficient(1:minMode),rightObj.Coefficient(end:-1:end - minMode + 1)); % Returns only the M^th coefficient of L*R
                else
                    L = leftObj.intlabpoly();
                    R = rightObj.intlabpoly();
                    LR = L*R;
                    switch truncationSize
                        case 'Fixed'
                            minMode = min([leftObj.Truncation,rightObj.Truncation]);
                            productObj = Scalar(LR, basis, minMode);
                        case 'Full'
                            productObj = Scalar(LR, basis);
                        otherwise
                            productObj = Scalar(LR, basis, truncationSize);
                    end
                end
            else
                switch truncationSize
                    case 'Fixed'
                        productSize = leftObj.Truncation;
                        productCoefficient = conv([zeros(productSize - 1, 1); leftObj.Coefficient],rightObj.Coefficient,'valid');
                        productObj = Scalar(productCoefficient, basis, productSize,Dimension);
                    case 'Recursion'
                        productObj = dot(leftObj.Coefficient,flip(rightObj.Coefficient)); % Returns only the M^th coefficient of L*R
                    case 'Full'
                        productObj = Scalar(conv(leftObj.Coefficient,rightObj.Coefficient), basis);
                    otherwise % specify an integer truncation
                        productCoefficient = conv(leftObj.Coefficient, rightObj.Coefficient);
                        productObj = Scalar(productCoefficient, basis, truncationSize);
                end
            end
        case 2 % ---------------------------------------- 2D SURFACES ----------------------------------
            if strcmp(leftObj.NumericalClass,'intval') || strcmp(rightObj.NumericalClass,'intval') % products for intval coefficients
                L = leftObj.intlabpoly();
                R = rightObj.intlabpoly();
                LR = L*R;
                switch truncationSize
                    case 'Fixed'
                        productObj = Scalar(LR,leftObj.Truncation, basis);
                    case 'Recursion'
                        minMode = min([leftObj.Truncation;rightObj.Truncation]);
                        mIdx = (LR.e(:,2) == minMode(1) - 1) & (LR.e(:,1) < minMode(2));
                        productObj = LR.c(mIdx)'; % Returns an intval vector, not a Scalar.
                    case 'Full'
                        productObj = Scalar(LR, basis, leftObj.Truncation + rightObj.Truncation - [1,1]);
                    otherwise
                        productObj = Scalar(LR, basis, truncationSize);
                end

            elseif isa(leftObj.Coefficient,'Scalar') && isa(rightObj.Coefficient,'Scalar')
                % compute coefficients of product
                if strcmp(truncationSize,'Full')
                    minIdx = 1;
                    maxIdx = leftObj.Truncation(1) + rightObj.Truncation(1) - 1; % full cauchy product (in 1st variable) has deg m + n - 1.
                    subTruncSize = 'Full';
                elseif strcmp(truncationSize,'Fixed') % Truncate Cauchy product to same size as factors in 1st variable and subTruncation in remainder.
                    minIdx = 1;
                    maxIdx = leftObj.Truncation(1);
                    subTruncSize = min([leftObj.Truncation(2:end),rightObj.Truncation(2:end)]);
                elseif strcmp(truncationSize,'Recursion')
                    minIdx = min([leftObj.Truncation(1),rightObj.Truncation(1)]);
                    maxIdx = minIdx;
                    subTruncSize = min([leftObj.Truncation(2:end),rightObj.Truncation(2:end)]);
                end

                % manually compute Cauchy product.
                for m = minIdx:maxIdx
                    lowerIdx = (1 + max([0,m - rightObj.Truncation(1)])):min([m,leftObj.Truncation(1)]);
                    upperIdx = m + 1 - lowerIdx;
                    leftContribution = leftObj.Coefficient(lowerIdx);
                    rightContribution = rightObj.Coefficient(upperIdx);
                    productCoefficient(m) = dot(leftContribution,rightContribution,subTruncSize);
                end

                if strcmp(truncationSize,'Fixed') || strcmp(truncationSize,'Full')
                    productObj = Scalar();
                    productObj.Basis = basis;
                    productObj.Dimension = leftObj.Dimension;
                    productObj.Coefficient = productCoefficient;
                    productTruncation = [productObj.Coefficient.Truncation];
                    productObj.Truncation = [maxIdx,max(productTruncation)];
                else
                    productObj = productCoefficient(m);
                end

            else % NumericalClass is double
                switch truncationSize
                    case 'Fixed'
                        productCoefficient = conv2([zeros(leftObj.Truncation-1),zeros(leftObj.Truncation-[1,0]);zeros(leftObj.Truncation-[0,1]),leftObj.Coefficient],rightObj.Coefficient,'valid');
                        % UPDATED FOR SCALAR CONSTRUCTOR productObj = Scalar(productCoefficient, basis, size(productCoefficient), Dimension);
                        productObj = Scalar(productCoefficient, basis, size(productCoefficient));

                    case 'FFT'
                        n = leftObj.Truncation + rightObj.Truncation - [1,1];
                        cauchyProduct = ifftn*(fftn(leftObj.Coefficient,n).*fftn(rightObj.Coefficient,n));
                        productObj = cauchyProduct(1:obj.Truncation(1),1:obj.Truncation(2));
                    case 'Recursion' % Object returned in this case is a double, not a Scalar.
                        productObj = conv2([zeros(size(leftObj.Coefficient)-[0,1]),leftObj.Coefficient],rightObj.Coefficient,'valid');
                    case 'Full'
                        productObj = Scalar(conv2(leftObj.Coefficient, rightObj.Coefficient), basis);
                    otherwise % truncationSize is a specfied size to truncate the product
                        if length(truncationSize) == leftObj.Dimension
                            productFull = mtimes(leftObj,rightObj,'Full');
                            productObj = Scalar(productFull.Coefficient, basis, truncationSize);
                        else
                            error('Multiplication of these Scalars is undefined')
                        end
                end
            end
        case 3 % ---------------------------------------- 3D SURFACES ----------------------------------
            switch truncationSize
                case 'Fixed'
                    modes = leftObj.Truncation;
                    padTo = leftObj.Truncation + rightObj.Truncation - [1,1,1];
                    fullCoefficient = ifftn(fftn(leftObj.Coefficient,padTo).*fftn(rightObj.Coefficient,padTo));
                    productObj = Scalar(fullCoefficient(1:modes(1), basis, 1:modes(2),1:modes(3)));
                case 'Recursion'
                    modes = leftObj.Truncation;
                    padTo = leftObj.Truncation + rightObj.Truncation - [1,1,1];
                    fullCoefficient = ifftn(fftn(leftObj.Coefficient,padTo).*fftn(rightObj.Coefficient,padTo));
                    productObj = Scalar(squeeze(fullCoefficient(modes(1), basis, 1:modes(2),1:modes(3))));
            end
    end
elseif isa(leftObj,'Scalar') % Scalar multiplication with Fscalar (right)
    basis = leftObj.Basis;
    productObj = Scalar(leftObj.Coefficient*rightObj, basis, leftObj.Truncation);
elseif isa(rightObj,'Scalar')% Scalar multiplication with Fscalar (left)
    basis = rightObj.Basis;

    try
        productObj = Scalar(leftObj*rightObj.Coefficient, basis, rightObj.Truncation);
    catch ME
        if strcmp(ME.message, 'Undefined function ''times'' for input arguments of type ''Scalar''.')
            disp('Left multiplication by intval is undefined. Use right multiplication only.')
        end
    end
end
end % mtimes

% Revision History:
%{
24-Jun-2017 - support for 1-dimensional Scalars
02-Jul-2017 - change multiplication to FFT based convolution
19-Aug-2017 - reverted back to non-FFT multiplication. Rigorous and numeric FFT needs to be impleneted more carefully before it can be trusted.
13-Aug-2018 - updated for Scalar class
%}
