function productObj = mtimes(leftObj,rightObj,varargin)
% Define multiplication of Scalars with Scalars, Fscalars (double or intval) and BAoperators (multiplication in usual matrix/vector sense). Current implementation only supports dimension-1 surfaces.

% Written by S. Kepley 03/2017
% Added support for 1D algebra 06/2017
% Changed intval products to use FFT 07/2017
% Reverted back to non-FFT multiplication. Rigorous FFT needs to be impleneted more carefully before using 08/2017

% ---------------------- INPUT ----------------------
% leftObj Scalar or BAoperator or double or intval: left factor
% rightObj Scalar or BAoperator or double or intval: right factor
% varargin{1} = truncation type (string): Use to allow fast multiplication via FFT % when full products aren't needed

% ---------------------- OUTPUT ----------------------
% productObj: Scalar corresponding to leftObj*rightObj


if isa(leftObj,'BAoperator') % BAoperator acts on Scalar (left only)
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

    % this is currently a catch all that needs to be fixed
    if isequal(leftObj.Truncation,[]) || isequal(rightObj.Truncation,[])
        disp('Scalar should never have empty value for Truncation')
        L = Scalar(leftObj.Coefficient);
        R = Scalar(rightObj.Coefficient);
        productObj = Scalar(mtimes(L,R,varargin{:}));
        productObj.Truncation = [];

    % one factor is a constant function
    elseif (leftObj.Dimension == 0 || rightObj.Dimension == 0)
        productObj = Scalar(leftObj.Coefficient*rightObj.Coefficient);

    % products require factors of the same dimension
    elseif ~isequal(leftObj.Dimension,rightObj.Dimension)
        error('Cauchy products for factors of different dimensions is not defined yet')

    % left/right surfaceDimension match. Set product surfaceDimension equal
    else
        surfaceDimension = leftObj.Dimension;
    end

    if (nargin > 2); % If both factors have double NumericalClass default behavior is to truncate the product to the same size as obj. Otherwise call with varargin.
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

    switch surfaceDimension
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
                            productObj = Scalar(LR,minMode);
                        case 'Full'
                            productObj = Scalar(LR);
                        otherwise
                            productObj = Scalar(LR,truncationSize);
                    end
                end
            else
                switch truncationSize
                    case 'Fixed'
                        productSize = leftObj.Truncation;
                        productCoefficient = conv([zeros(1,productSize - 1),leftObj.Coefficient],rightObj.Coefficient,'valid');
                        productObj = Scalar(productCoefficient,productSize,surfaceDimension);
                    case 'Recursion'
                        productObj = dot(leftObj.Coefficient,flip(rightObj.Coefficient)); % Returns only the M^th coefficient of L*R
                    case 'Full'
                        productObj = Scalar(conv(leftObj.Coefficient,rightObj.Coefficient));
                    otherwise % specify an integer truncation
                        productCoefficient = conv(leftObj.Coefficient,rightObj.Coefficient);
                        productObj = Scalar(productCoefficient,truncationSize);
                end
            end
        case 2 % ---------------------------------------- 2D SURFACES ----------------------------------
            if strcmp(leftObj.NumericalClass,'intval') || strcmp(rightObj.NumericalClass,'intval') % products for intval coefficients
                L = leftObj.intlabpoly();
                R = rightObj.intlabpoly();
                LR = L*R;
                switch truncationSize
                    case 'Fixed'
                        productObj = Scalar(LR,leftObj.Truncation);
                    case 'Recursion'
                        minMode = min([leftObj.Truncation;rightObj.Truncation]);
                        mIdx = (LR.e(:,2) == minMode(1) - 1) & (LR.e(:,1) < minMode(2));
                        productObj = LR.c(mIdx)'; % Returns an intval vector, not a Scalar.
                    case 'Full'
                        productObj = Scalar(LR,leftObj.Truncation + rightObj.Truncation - [1,1]);
                    otherwise
                        productObj = Scalar(LR,truncationSize);
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
                        productObj = Scalar(productCoefficient,size(productCoefficient),surfaceDimension);
                    case 'FFT'
                        n = leftObj.Truncation + rightObj.Truncation - [1,1];
                        cauchyProduct = ifftn*(fftn(leftObj.Coefficient,n).*fftn(rightObj.Coefficient,n));
                        productObj = cauchyProduct(1:obj.Truncation(1),1:obj.Truncation(2));
                    case 'Recursion' % Object returned in this case is a double, not a Scalar.
                        productObj = conv2([zeros(size(leftObj.Coefficient)-[0,1]),leftObj.Coefficient],rightObj.Coefficient,'valid');
                    case 'Full'
                        productObj = Scalar(conv2(leftObj.Coefficient,rightObj.Coefficient));
                    otherwise % truncationSize is a specfied size to truncate the product
                        if length(truncationSize) == leftObj.Dimension;
                            productFull = mtimes(leftObj,rightObj,'Full');
                            productObj = Scalar(productFull.Coefficient,truncationSize);
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
                    productObj = Scalar(fullCoefficient(1:modes(1),1:modes(2),1:modes(3)));
                case 'Recursion'
                    modes = leftObj.Truncation;
                    padTo = leftObj.Truncation + rightObj.Truncation - [1,1,1];
                    fullCoefficient = ifftn(fftn(leftObj.Coefficient,padTo).*fftn(rightObj.Coefficient,padTo));
                    productObj = Scalar(squeeze(fullCoefficient(modes(1),1:modes(2),1:modes(3))));
            end
    end
elseif isa(leftObj,'Scalar') % Scalar multiplication with Fscalar (right)
    productObj = Scalar(leftObj.Coefficient*rightObj,leftObj.Truncation);
elseif isa(rightObj,'Scalar')% Scalar multiplication with Fscalar (left)
    try
        productObj = Scalar(leftObj*rightObj.Coefficient,rightObj.Truncation);
    catch ME
        if strcmp(ME.message, 'Undefined function ''times'' for input arguments of type ''Scalar''.')
            disp('Left multiplication by intval is undefined. Use right multiplication only.')
        end
    end
end
end
