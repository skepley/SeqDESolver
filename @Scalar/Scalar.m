classdef Scalar
%SCALAR - SeqDE class for representing analytic scalar functions.
%
%   The Scalar class is a finite approximation for an analytic scalar of the form, f: D --> R^n where D is
%   a d-dimensional polydisc in C^d. Scalar representations are specified as double or intval coefficients
%   with respect to a Taylor, Fourier, or Chebyshev series expansion.
%
%   SCALAR constructor syntax:
%       ScalarObj = Scalar()
%
%   SCALAR properties:
%       Coefficient - array of coefficients
%       Dimension - dimension of the polydisc
%       Truncation - vector of truncation degrees for each dimension
%       Weight - vector of weights for sequence space (default is [1,...,1])
%       NumericalClass - string specifying double, interval, or Scalar coefficients
%
%   SCALAR methods:
%       intval - convert to interval coefficients
%
%   Examples:
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
%   Subclasses: none
%   Superclasses: none
%   Other classes required: none
%   Other m-files required: INTLAB toolbox
%   MAT-files required: none
%
%   See also: @BAscalar (previous version of this class), @Chart

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Jun-2016; Last revision: 08-Aug-2018
%
%   ToDo:
%   Merge and simplify evaluation, fixtime, derivative, and truncation methods/function calls.
%   Rigorous evaluation of sin/cos/exp, inverses, and fractional powers.
%   Finish vectorization of all methods.
%   Finish high precision inteval FFT and reimplement FFT-based fast convolution.
%   Add polynomial evaluation over the Banach algebra.
%   Fix shift method to use new subscript reference.


    %% -------------------- Properties --------------------
    properties
        Basis; % Taylor, Fourier, Chebyshev: string
        Coefficient; % N1-by-N2-by-...-Nd array: double or intval.
        Dimension; % d: integer
        NumericalClass; % double, intval, Scalar: string
        Truncation; % [N1,...,Nd]: integer
    end



    properties(Hidden = 1)
        Weight = 'ones'; % [nu_1,...,nu_d]: positive double

    end


    %% -------------------- Methods --------------------
    methods
        function obj = Scalar(coefficient, varargin)
        %SCALAR - class constructor

            if nargin > 0

                if nargin > 2;
                    obj.Dimension = varargin{2}; % specify dimension as vararg.
                end

                if nargin > 3; % support for weighted ell^1 vectors with arbitrary weights.
                    obj.Weight = varargin{3};
                end

                %% ----------------------------------- Coefficient specifed as cell array  -----------------------------------
                if isa(Coefficient,'cell') % cell array of coefficients ---> vector of Scalars
                    if nargin ==1; % specify a vector of Scalars with fixed dimension specified by varargin{1}.
                        objLength = length(Coefficient);
                        obj(objLength) = Scalar(Coefficient{objLength},varargin{:});
                        for j = 1:objLength-1
                            obj(j) = Scalar(Coefficient{j},varargin{:});
                        end
                        obj = reshape(obj,size(Coefficient));

                    elseif isa(varargin{1},'cell') % specify a single Scalar with Scalar coefficients specified by the cells of Coefficient.
                        % In this case obj.Truncation is a single integer and obj.Coefficient is a vector of Scalars with non-uniform number of modes.
                        obj.NumericalClass = 'Scalar';
                        switch numel(varargin{1})

                            case 1
                                obj.NumericalClass = 'Scalar';
                                obj.Dimension = varargin{2}; % must specify surface Dimension explicitly.
                                obj.Truncation = [varargin{1},Inf(1,obj.Dimension - 1)];
                                surfaceCoefficient(obj.Truncation(1)) = Scalar(); % initiate coefficient vector.
                                for j = 1:length(Coefficient)
                                    surfaceCoefficient(j) = Scalar(Coefficient{j},Inf(1,obj.Dimension-1));
                                end

                            case 2 % {M,[N1,N2,...]} for fixed [N1,N2,...] modes for (d-1) dimensional coefficients
                                subTruncation = varargin{1}{2};
                                obj.Truncation = [varargin{1},subTruncation];
                                obj.Dimension = length(obj.Truncation);
                                obj.Coefficient(obj.Truncation(1)) = Scalar(); % initiate vector of (d-1)-dimensional Scalar
                                for j = 1:length(Coefficient)
                                    obj.Coefficient(j) = Scalar(Coefficient{j},subTruncation);
                                end
                        end

                    end

                %% ----------------------------------- Coefficient specifed as intval polynomial  -----------------------------------
                elseif isa(Coefficient,'polynom')
                    if nargin > 1
                        obj.Truncation = varargin{1};
                    else
                        obj.Truncation = Coefficient.e + ones(1,length(Coefficient.e));
                    end
                    obj.Dimension = length(obj.Truncation); % length of Degree equals dimension of surface

                    switch obj.Dimension
                        case 1
                            intvalCoefficient = flip(Coefficient.c);
                            obj.Coefficient = intvalCoefficient(1:min(end,obj.Truncation));
                        case 2
                            deg = max(Coefficient.e);
                            intvalCoefficient = reshape(flip(flip(Coefficient.c,1),2),1 + deg(2),[]); % full product coefficients
                            obj.Coefficient = intvalCoefficient(1:min(end,obj.Truncation(1)),1:min(end,obj.Truncation(2)));
                        otherwise
                            error('Not yet implemented')
                    end
                    obj.NumericalClass = class(obj.Coefficient);

                %% ----------------------------------- Coefficient specifed as Scalar  -----------------------------------
                elseif isa(Coefficient,'Scalar') % copy to a new Scalar
                    if nargin == 1;
                        obj.Coefficient = Coefficient.Coefficient;
                        obj.Truncation = Coefficient.Truncation;
                        obj.NumericalClass = Coefficient.NumericalClass;
                        obj.Dimension = Coefficient.Dimension;
                    elseif nargin ==2
                        obj = Scalar(Coefficient.Coefficient,varargin{1});
                    end

                %% ----------------------------------- Coefficient specifed as double or intval array  -----------------------------------
                else
                    switch nargin
                        case 1 % input is Coefficient of correct size
                            obj.Coefficient = Coefficient;
                            obj.NumericalClass = class(obj.Coefficient);
                            dims = size(Coefficient);
                            trueDims = dims(dims > 1); % remove singleton dimensions
                            obj.Dimension = length(trueDims);
                            if obj.Dimension > 0
                                obj.Truncation = trueDims;
                            else
                                obj.Truncation = 1;
                            end

                        otherwise % input is Coefficient and truncation
                            if isa(varargin{1},'cell')
                                if numel(varargin{1}) == 1 % {M} specify only modes in time. coefficient type is Scalar with flexible modes for coefficients
                                    obj.NumericalClass = 'Scalar';
                                    obj.Dimension = varargin{2}; % must specify surface Dimension explicitly.
                                    obj.Truncation = [varargin{1},Inf(1,obj.Dimension - 1)];

                                    surfaceCoefficient(obj.Truncation(1)) = Scalar(); % initiate coefficient vector.
                                    switch obj.Dimension
                                        case 2
                                            for j = 1:size(Coefficient,1)
                                                surfaceCoefficient(j) = Scalar(Coefficient(j,:));
                                            end
                                        case 3
                                            for j = 1:size(Coefficient,1)
                                                surfaceCoefficient(j) = Scalar(Coefficient(j,:,:));
                                            end
                                        otherwise
                                            error('Scalar coefficients supported for dimension 2 or 3 only')
                                    end
                                else % coefficient type is Scalar with uniform modes for coefficients
                                    obj.NumericalClass = 'Scalar';
                                    subTruncation = varargin{1}{2}; % [N1,N2,...] degree of d-1 dimensional coefficients
                                    obj.Truncation = [varargin{1}{1},subTruncation]; % [M,N1,N2,...]
                                    obj.Dimension = length(subTruncation) + 1;
                                    surfaceCoefficient(obj.Truncation(1)) = Scalar(); % initiate vector of (d-1)-dimensional Scalar
                                    switch obj.Dimension
                                        case 2
                                            for j = 1:size(Coefficient,1)
                                                surfaceCoefficient(j) = Scalar(Coefficient(j,:),subTruncation);
                                            end
                                        case 3
                                            for j = 1:size(Coefficient,1)
                                                surfaceCoefficient(j) = Scalar(Coefficient(j,:,:),subTruncation);
                                            end
                                        otherwise
                                            error('Scalar coefficients supported for dimension 2 or 3 only')
                                    end
                                    obj.Coefficient = surfaceCoefficient;
                                end
                            else % intval or double Coefficient with modes specified as double
                                obj.Truncation = varargin{1};
                                obj.Dimension = length(obj.Truncation); % length of Degree equals dimension of surface
                                if isa(Coefficient,'double') || isa(Coefficient,'intval') % coefficients given as double or intval array
                                    switch obj.Dimension
                                        case 0
                                            obj.Coefficient = Coefficient;
                                        case 1
                                            obj.Coefficient = Coefficient(1:min(end,obj.Truncation));
                                        case 2
                                            obj.Coefficient = Coefficient(1:min(end,obj.Truncation(1)),1:min(end,obj.Truncation(2)));
                                        case 3
                                            if obj.Truncation(1) == 1
                                                obj.Coefficient = Coefficient(1:min(end,obj.Truncation(2)),1:min(end,obj.Truncation(3)));
                                            else
                                                obj.Coefficient = Coefficient(1:min(end,obj.Truncation(1)),1:min(end,obj.Truncation(2)),1:min(end,obj.Truncation(3)));
                                            end
                                        otherwise
                                            error('Not yet implemented')
                                    end
                                    obj.NumericalClass = class(obj.Coefficient);
                                elseif isa(Coefficient,'Scalar')
                                    if length(Coefficient) == 1
                                        obj = Scalar(Coefficient.Coefficient,obj.Truncation);
                                    else
                                        obj = Coefficient;
                                    end
                                end
                                obj = padcoefficient(obj)
                            end
                    end
                end
            end
        end % end class constructor



        function paddedObj = padcoefficient(obj)
        %PADCOEFFICIENT - pads Scalar coefficient with zeros to achieve truncation size consistent with obj.Truncation

            if length(obj) > 1
                return
            elseif isequal(obj.Truncation, size(obj.Coefficient)) % already padded
                return
            end

            switch obj.Dimension
                case 0
                    paddedObj = obj;
                    return

                case 1
                    if isequal(length(obj.coef), obj.Truncation)
                        return
                    elseif strcmp(obj.NumericalClass,'double')
                        coef = zeros(1,obj.Truncation);
                    elseif strcmp(obj.NumericalClass,'intval')
                        coef = midrad(zeros(1,obj.Truncation),0);
                    end
                    coef(1:length(obj.Coefficient)) = obj.Coefficient;

                case 2
                    if isequal(size(obj.coef),obj.Truncation)
                        return
                    elseif strcmp(obj.NumericalClass,'double')
                        coef = zeros(obj.Truncation);
                    elseif strcmp(obj.NumericalClass,'intval')
                        coef = midrad(zeros(obj.Truncation),0);
                    end
                    coefSize = size(obj.Coefficient);
                    coef(1:coefSize(1),1:coefSize(2)) = obj.Coefficient;

                case 3
                    if isequal(size(obj.coef),obj.Truncation)
                        return
                    elseif strcmp(obj.NumericalClass,'double')
                        coef = zeros(obj.Truncation);
                    elseif strcmp(obj.NumericalClass,'intval')
                        coef = midrad(zeros(obj.Truncation),0);
                    end
                    coefSize = size(obj.Coefficient);
                    coef(1:coefSize(1),1:coefSize(2),1:coefSize(3)) = obj.Coefficient;

                otherwise
                    error('padcoefficient not implemented for this dimension')
            end
            obj.Coefficient = coef;
            paddedObj = obj;
        end



        %% -------------------- METHOD REPAIR SHOP --------------------
        % The following methods need to be updated and thoroughly checked before using.

        function int_ds = intds(obj,varargin)
            % evaluate definite or indefinite integral with respect to spatial variable
            if obj.Dimension == 2
                C = repmat(1./(1:obj.Truncation(2)),obj.Truncation(1),1);
                int_ds = Scalar([zeros(obj.Truncation(1),1),C.*obj.Coefficient]); % indefinite integral
                if (nargin > 1) % compute integral obj ds on [a,b]
                    bounds = varargin{1};
                    int_ds = int_ds.fix_space(bounds(2)) - int_ds.fix_space(bounds(1)); % evaluation definite integral
                end
            else
                error('ds not implemented for this dimension')
            end
        end


        % THIS NEEDS TO BE IMPLEMENTED USING INTLAB POLYEVAL AND MERGED INTO THE EVAL METHOD.
        function evalObj = intvalEval(obj,s,t)
            % s is a vector of intervals, t is a single scalar.
            X = obj.fixtime(t);
            F = polynom(flip(X.Coefficient),'s');
            evalObj = polyval(F,s);
        end

        % THIS NEEDS TO BE MERGED INTO THE EVAL METHOD.
        function evalObj = fixSpace(obj,s)
            % collapse Taylor series onto fixed space evaluation, s
            space_vector = bsxfun(@power,s,(0:size(obj.Coefficient,2)-1)');
            evalObj = obj.Coefficient*space_vector;
        end





        %% -------------------- METHOD GRAVEYARD --------------------
        % The following methods should not be called but are retained for a while to make sure they don't break any old code.

        function newObj = homog(obj)
           % convert Scalar coefs to double or intval coefs

            warning('homog is deprecated and Scalar NumercalClass should never be used')

           arrayCoefficient = [obj.Coefficient.Coefficient];
           newCoefficient = reshape(arrayCoefficient, obj.Truncation(1), []);
           newObj = Scalar(newCoefficient', obj.Truncation, obj.Dimension);
        end



        function truncate(obj,modes)
            % truncate Scalar to specified number of modes

            warning('truncate method is deprecated. It should be carefully checked before trusting results')
            if ~isequal(obj.Dimension,2)
                error('not implemented for dimension other than 2 yet')
            end

            if isequal(length(obj),1) % a single Scalar
                newTruncation = [min([obj.Truncation(1),modes(1)]),min([obj.Truncation(2),modes(2)])];
                obj.Coefficient = obj.Coefficient(1:newTruncation(1),1:newTruncation(2));
                obj.Coefficient = obj.Coefficient(1:newTruncation(1),1:newTruncation(2));
                obj.Coefficient = obj.Coefficient(1:newTruncation(1),1:newTruncation(2));
                obj.Coefficient = obj.Coefficient(1:newTruncation(1),1:newTruncation(2));
                obj.Truncation = newTruncation;
            else
                for j = 1:length(obj)
                   truncate(obj(j),modes)
                end
            end
        end


        function columnObj = col(obj)
            % returns Coefficients of BAsclar as a column vector (double) under the canonical isomorphism

            warning('col is deprecated and should be replaced by the exponent method')
            columnObj = reshape(obj.Coefficient,[],1);
        end


        function coefArray = coef(obj)
            % returns a double or intval array of the coefficients of obj

            warning('coef is deprecated and should be replaced by the exponent method')
            if length(obj) == 1
                coefArray = squeeze(obj.Coefficient);
            else
                cellCoefficient = arrayfun(@(j)obj(j).Coefficient,1:length(obj),'UniformOutput',false);
                % coefArray = cell2mat(cellCoefficient');
                coefArray = cellCoefficient;
            end
        end

        function dObj_ds = ds(obj)
            % compute spatial derivative

            warning('ds is deprecated. Use dt method instead.')
            switch obj.Dimension
                case 1
                    C = 1:obj.Truncation-1;
                    dObj_ds = Scalar(C.*obj.Coefficient(2:end));
                case 2
                    C = repmat((1:obj.Truncation(2)-1),obj.Truncation(1),1);
                    dObj_ds = Scalar(C.*obj.Coefficient(:,2:obj.Truncation(2)));
                otherwise
                    error('ds not implemented for this dimension')
            end
        end

        function evalObj = grideval(obj,varargin)

            warning('grideval is deprecated and should be replaced by the eval method')
            switch nargin
                case 3
                    s = varargin{1};
                    t = varargin{2};
                    % evaluation is in meshgrid format. This should be changed to ndgrid and incorporated into the eval method.

                    flipCoefficients = fliplr(flip(obj.Coefficient)); %switch Coefficient to descending powers
                    evalSpatial = nan(length(s),obj.Truncation(1));
                    for j = 1:obj.Truncation(1)
                        evalSpatial(:,j) = polyval(flipCoefficients(j,:),s);
                    end
                    evalObj = nan(length(t),length(s));
                    for k = 1:length(s)
                        evalObj(:,k) = polyval(evalSpatial(k,:),t);
                    end

                case 4 % 2d evaluation is in ndgrid format
                    % TEST EVALUATION
                    % modes = [3,3,3]
                    % A = zeros(modes);
                    % A(1,2,1) = 1;
                    % A(1,1,2) = 1;
                    % A(2,1,1) = 1;
                    % A corresponds to P(x,y,z) = x + y + z
                    %
                    % a = Scalar(A);
                    % s1 = [1,-2,3]
                    % s2 = [0,1]
                    % t = [2,3,4,5]
                    % checkme = a.gridEval(s1,s2,t)

                    s1 = varargin{1};
                    s2 = varargin{2};
                    t = varargin{3};

                    evalDims = [length(s1),length(s2),length(t)];
                    evalObj = nan(evalDims);
                    coef = zeros(length(s2),length(s1),obj.Truncation(1));
                    for j = 1:obj.Truncation(1)
                        pj = Scalar(squeeze(obj.Coefficient(j,:,:)));
                        coef(:,:,j) = pj.gridEval(s1,s2);
                    end

                    for k = 1:evalDims(1)
                        for l = 1:evalDims(2)
                            evalObj(k,l,:) = polyval(flip(squeeze(coef(l,k,:))),t);
                        end
                    end
            end
        end

    end % end methods

    %% STATIC METHODS
    methods(Static)
        function zarray = zeros(varargin)
            if nargin ==0
                zarray = Scalar(0,[5,5]);
            else
                zarray = repmat(Scalar(0,[5,5]),varargin{:});
            end
        end
    end % end static methods
end % end classdef

% Revision History:
%{
11-Jul-2017 - Support for interval and Scalar coefficients added.
15-Aug-2017 - Reverted FFT based convolution to a classical algorithm. FFT is faster but numerically unstable, especially for intval
    coefficients.
08-Aug-2018 - Class completely overhauled and renamed from BAscalar to SCalar. Numerous improvements including:
   full class code refactorization and organization of class folder
   class inheritance changed from handle to value
   support for additional bases (Fourier or Chebyshev)
   Intlab polynom dependency removed
   support for higher dimensions
   streamlined truncation methods
   subscript reference methods
   numerous efficiency gains, bug fixes, and improvements in data type consistency
%}


