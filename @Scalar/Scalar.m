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
    %       coefficient - array of coefficients
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
    %   See also: @Scalar (previous version of this class), @Chart

    %   Author: Shane Kepley
    %   email: shane.kepley@rutgers.edu
    %   Date: 08-Jun-2016; Last revision: 18-Aug-2018
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
        function obj = Scalar(coefficient, basis, varargin)
            %SCALAR - class constructor

            if nargin > 0
                if nargin > 3
                    obj.Dimension = varargin{2}; % specify dimension as vararg.
                end

                if nargin > 4 % support for weighted ell^1 vectors with arbitrary weights.
                    obj.Weight = varargin{3};
                end

                obj.Basis = basis;
                
                if sum(size(coefficient) > 1) == 1 % reshape any 1-d arrays into a column vector
                    coefficient = reshape(coefficient, [],1); %
                end

                %% ----------------------------------- coefficient specifed as cell array  -----------------------------------
                if isa(coefficient,'cell') % cell array of coefficients ---> vector of Scalars
                    if nargin ==1 % specify a vector of Scalars with fixed dimension specified by varargin{1}.
                        objLength = length(coefficient);
                        obj(objLength) = Scalar(coefficient{objLength}, basis, varargin{:});
                        for j = 1:objLength-1
                            obj(j) = Scalar(coefficient{j}, basis, varargin{:});
                        end
                        obj = reshape(obj, size(coefficient));

                    elseif isa(varargin{1},'cell') % specify a single Scalar with Scalar coefficients specified by the cells of coefficient.
                        % In this case obj.Truncation is a single integer and obj.Coefficient is a vector of Scalars with non-uniform number of modes.
                        obj.NumericalClass = 'Scalar';
                        switch numel(varargin{1})

                            case 1
                                error('Scalar - this should not be trusted')
                                obj.NumericalClass = 'Scalar';
                                obj.Dimension = varargin{2}; % must specify surface Dimension explicitly.
                                obj.Truncation = [varargin{1}, Inf(1,obj.Dimension - 1)];
                                surfacecoefficient(obj.Truncation(1)) = Scalar(); % initiate coefficient vector.
                                for j = 1:length(coefficient)
                                    surfacecoefficient(j) = Scalar(coefficient{j}, basis, Inf(1,obj.Dimension-1));
                                end

                            case 2 % {M,[N1,N2,...]} for fixed [N1,N2,...] modes for (d-1) dimensional coefficients
                                subTruncation = varargin{1}{2};
                                obj.Truncation = [varargin{1}, subTruncation];
                                obj.Dimension = length(obj.Truncation);
                                obj.Coefficient(obj.Truncation(1)) = Scalar(); % initiate vector of (d-1)-dimensional Scalar
                                for j = 1:length(coefficient)
                                    obj.Coefficient(j) = Scalar(coefficient{j}, basis, subTruncation);
                                end
                        end

                    end

                    %% ----------------------------------- coefficient specifed as intval polynomial  -----------------------------------
                elseif isa(coefficient,'polynom')
                    if nargin > 1
                        obj.Truncation = varargin{1};
                    else
                        obj.Truncation = coefficient.e + ones(1,length(coefficient.e));
                    end
                    obj.Dimension = length(obj.Truncation); % length of Degree equals dimension of surface

                    switch obj.Dimension
                        case 1
                            intvalcoefficient = flip(coefficient.c);
                            obj.Coefficient = intvalcoefficient(1:min(end,obj.Truncation));
                        case 2
                            deg = max(coefficient.e);
                            intvalcoefficient = reshape(flip(flip(coefficient.c,1),2),1 + deg(2),[]); % full product coefficients
                            obj.Coefficient = intvalcoefficient(1:min(end,obj.Truncation(1)),1:min(end,obj.Truncation(2)));
                        otherwise
                            error('Not yet implemented')
                    end
                    obj.NumericalClass = class(obj.Coefficient);

                    %% ----------------------------------- coefficient specifed as Scalar  -----------------------------------
                elseif isa(coefficient,'Scalar') % copy to a new Scalar
                    if nargin == 2
                        obj.Coefficient = coefficient.coefficient;
                        obj.Truncation = coefficient.Truncation;
                        obj.NumericalClass = coefficient.NumericalClass;
                        obj.Dimension = coefficient.Dimension;
                    elseif nargin == 3
                        obj = Scalar(coefficient.coefficient, coefficient.Basis, varargin{1});
                    end

                    %% ----------------------------------- coefficient specifed as double or intval array  -----------------------------------
                else
                    switch nargin
                        case 2 % input is coefficient of correct size
                            obj.Coefficient = coefficient;
                            obj.NumericalClass = class(obj.Coefficient);
                            dims = size(coefficient);
                            trueDims = dims(dims > 1); % remove singleton dimensions
                            obj.Dimension = length(trueDims);
                            if obj.Dimension > 0
                                obj.Truncation = trueDims;
                            else
                                obj.Truncation = 1;
                            end

                        otherwise % input is coefficient and truncation
                            if isa(varargin{1},'cell')
                                if numel(varargin{1}) == 1 % {M} specify only modes in time. coefficient type is Scalar with flexible modes for coefficients
                                    obj.NumericalClass = 'Scalar';
                                    obj.Dimension = varargin{2}; % must specify surface Dimension explicitly.
                                    obj.Truncation = [varargin{1},Inf(1,obj.Dimension - 1)];

                                    surfacecoefficient(obj.Truncation(1)) = Scalar(); % initiate coefficient vector.
                                    switch obj.Dimension
                                        case 2
                                            for j = 1:size(coefficient,1)
                                                surfacecoefficient(j) = Scalar(coefficient(j,:), basis);
                                            end
                                        case 3
                                            for j = 1:size(coefficient,1)
                                                surfacecoefficient(j) = Scalar(coefficient(j,:,:), basis);
                                            end
                                        otherwise
                                            error('Scalar coefficients supported for dimension 2 or 3 only')
                                    end
                                else % coefficient type is Scalar with uniform modes for coefficients
                                    obj.NumericalClass = 'Scalar';
                                    subTruncation = varargin{1}{2}; % [N1,N2,...] degree of d-1 dimensional coefficients
                                    obj.Truncation = [varargin{1}{1},subTruncation]; % [M,N1,N2,...]
                                    obj.Dimension = length(subTruncation) + 1;
                                    surfacecoefficient(obj.Truncation(1)) = Scalar(); % initiate vector of (d-1)-dimensional Scalar
                                    switch obj.Dimension
                                        case 2
                                            for j = 1:size(coefficient,1)
                                                surfacecoefficient(j) = Scalar(coefficient(j,:), basis, subTruncation);
                                            end
                                        case 3
                                            for j = 1:size(coefficient,1)
                                                surfacecoefficient(j) = Scalar(coefficient(j,:,:), basis, subTruncation);
                                            end
                                        otherwise
                                            error('Scalar coefficients supported for dimension 2 or 3 only')
                                    end
                                    obj.Coefficient = surfacecoefficient;
                                end
                            else % intval or double coefficient with modes specified as double
                                obj.Truncation = varargin{1};
                                obj.Dimension = length(obj.Truncation); % length of Degree equals dimension of surface
                                if isa(coefficient,'double') || isa(coefficient,'intval') % coefficients given as double or intval array
                                    switch obj.Dimension
                                        case 0
                                            obj.Coefficient = coefficient;
                                        case 1
                                            obj.Coefficient = coefficient(1:min(end,obj.Truncation));
                                        case 2
                                            obj.Coefficient = coefficient(1:min(end,obj.Truncation(1)),1:min(end,obj.Truncation(2)));
                                        case 3
                                            if obj.Truncation(1) == 1
                                                obj.Coefficient = coefficient(1:min(end,obj.Truncation(2)),1:min(end,obj.Truncation(3)));
                                            else
                                                obj.Coefficient = coefficient(1:min(end,obj.Truncation(1)),1:min(end,obj.Truncation(2)),1:min(end,obj.Truncation(3)));
                                            end
                                        otherwise
                                            error('Not yet implemented')
                                    end
                                    obj.NumericalClass = class(obj.Coefficient);
                                elseif isa(coefficient,'Scalar')
                                    if length(coefficient) == 1
                                        obj = Scalar(coefficient.coefficient, basis, obj.Truncation);
                                    else
                                        obj = coefficient;
                                    end
                                end % if isa(coefficient, *)
                                obj = padcoefficient(obj);
                            end % if isa(varargin{1}, *)
                    end % switch nargin
                end % isa(*, coefficient)
            end % if
        end %  class constructor



        function paddedObj = padcoefficient(obj)
            %PADCOEFFICIENT - pads Scalar coefficient with zeros to achieve truncation size consistent with obj.Truncation

            if length(obj) > 1
                return
            elseif isequal(obj.Truncation, size(obj.Coefficient)) % already padded
                return
            end

            switch obj.Dimension
                case 0 % constant
                    paddedObj = obj;
                    return

                case 1 % 1-d surface
                    if isequal(length(obj.Coefficient), obj.Truncation)
                        paddedObj = obj;
                        return
                    elseif strcmp(obj.NumericalClass,'double')
                        coefficient = zeros(obj.Truncation, 1);
                    elseif strcmp(obj.NumericalClass,'intval')
                        coefficient = intval(zeros(obj.Truncation, 1),0);
                    end
                    coefficient(1:length(obj.Coefficient)) = obj.Coefficient;

                case 2 % 2-d surface
                    if isequal(size(obj.Coefficient),obj.Truncation)
                        return
                    elseif strcmp(obj.NumericalClass,'double')
                        coefficient = zeros(obj.Truncation);
                    elseif strcmp(obj.NumericalClass,'intval')
                        coefficient = midrad(zeros(obj.Truncation),0);
                    end
                    coefSize = size(obj.Coefficient);
                    coefficient(1:coefSize(1),1:coefSize(2)) = obj.Coefficient;

                case 3 % 3-d surface
                    if isequal(size(obj.Coefficient),obj.Truncation)
                        return
                    elseif strcmp(obj.NumericalClass,'double')
                        coefficient = zeros(obj.Truncation);
                    elseif strcmp(obj.NumericalClass,'intval')
                        coefficient = midrad(zeros(obj.Truncation),0);
                    end
                    coefSize = size(obj.Coefficient);
                    coefficient(1:coefSize(1),1:coefSize(2),1:coefSize(3)) = obj.Coefficient;

                otherwise
                    error('padcoefficient not implemented for this dimension')
            end
            obj.Coefficient = coefficient;
            paddedObj = obj;
        end



        %% -------------------- METHOD REPAIR SHOP --------------------
        % The following methods need to be updated and thoroughly checked before using.

        function int_ds = intds(obj, varargin)
            % evaluate definite or indefinite integral with respect to spatial variable
            if obj.Dimension == 2
                C = repmat(1./(1:obj.Truncation(2)),obj.Truncation(1),1);
                int_ds = Scalar([zeros(obj.Truncation(1),1), basis, C.*obj.Coefficient]); % indefinite integral
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
            F = polynom(flip(X.coefficient),'s');
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

            arraycoefficient = [obj.Coefficient.coefficient];
            newcoefficient = reshape(arraycoefficient, obj.Truncation(1), []);
            newObj = Scalar(newcoefficient', obj.Basis,  obj.Truncation, obj.Dimension);
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
            % returns coefficients of BAsclar as a column vector (double) under the canonical isomorphism

            warning('col is deprecated and should be replaced by the exponent method')
            columnObj = reshape(obj.Coefficient,[],1);
        end


        function dObj_ds = ds(obj)
            % compute spatial derivative

            warning('ds is deprecated. Use dt method instead.')
            switch obj.Dimension
                case 1
                    C = 1:obj.Truncation-1;
                    dObj_ds = Scalar(C.*obj.Coefficient(2:end), obj.Basis);
                case 2
                    C = repmat((1:obj.Truncation(2)-1),obj.Truncation(1),1);
                    dObj_ds = Scalar(C.*obj.Coefficient(:,2:obj.Truncation(2)), obj.Basis);
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

                    flipcoefficients = fliplr(flip(obj.Coefficient)); %switch coefficient to descending powers
                    evalSpatial = nan(length(s),obj.Truncation(1));
                    for j = 1:obj.Truncation(1)
                        evalSpatial(:,j) = polyval(flipcoefficients(j,:),s);
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
                    coefficient = zeros(length(s2),length(s1),obj.Truncation(1));
                    for j = 1:obj.Truncation(1)
                        pj = Scalar(squeeze(obj.Coefficient(j,:,:)), obj.Basis);
                        coefficient(:,:,j) = pj.gridEval(s1,s2);
                    end

                    for k = 1:evalDims(1)
                        for l = 1:evalDims(2)
                            evalObj(k,l,:) = polyval(flip(squeeze(coefficient(l,k,:))),t);
                        end
                    end
            end
        end

    end %  methods

    %% STATIC METHODS
    methods(Static)
        function zarray = zeros(varargin)
            if nargin ==0
                zarray = Scalar(0,[5,5]);
            else
                zarray = repmat(Scalar(0,[5,5]),varargin{:});
            end
        end
    end %  static methods
end %  classdef

% Revision History:
%{
11-Jul-2017 - Support for interval and Scalar coefficients added.
15-Aug-2017 - Reverted FFT based convolution to a classical algorithm. FFT is faster but numerically unstable, especially for intval
    coefficients.
08-Aug-2018 - Class completely overhauled and renamed from Scalar to SCalar. Numerous improvements including:
   full class code refactorization and organization of class folder
   class inheritance changed from handle to value
   support for additional bases (Fourier or Chebyshev)
   Intlab polynom dependency removed
   support for higher dimensions
   streamlined truncation methods
   subscript reference methods
   numerous efficiency gains, bug fixes, and improvements in data type consistency
18-Aug-2018 - 1-dimensional Scalars changed to column vector coefficients to speed up operations
%}


