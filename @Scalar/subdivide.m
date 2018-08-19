function varargout = subdivide(obj, parmRange, varargin)
%SUBDIVIDE - rigorously subdivide a Scalar parameterization into a piecewise Scalar parameterization.
%
%   Syntax:
%       subScalar = SUBDIVIDE(obj, parmRange) is a vector of Scalars of length nDivs which is the number of rows of parmRange.
%           Each row of parmRange specifies a subset of D^d on which to parameterize each subScalar.
%
%       subScalar = SUBDIVIDE(obj, parmRange, t0) assumes that obj is a parameterization in time with respect to the first dimension. obj is first
%           evaluated at tau = t0 and the resulting (d-1)-dimensional Scalar is subdivided into subScalars.
%           a piecewise parameterizationsubdivided
%
%       [subScalar, subDomain] = SUBDIVIDE(obj, parmRange) subDomain is a nDivs-by-(2d) array of material coordinates with respect to [-1,1]^d. These coordinates
%           are consecutively inherited by successive calls of subdivide.
%
%   Inputs:
%       obj - Scalar
%       parmRange - nDivs-by-(2d) array
%       t0 - double
%
%   Outputs:
%       newParm - vector of Scalars
%       newMTC - double array
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jan-2017; Last revision: 13-Aug-2018

% parse input
p = inputParser;
addRequired(p,'obj')
addRequired(p,'parmRange')
addOptional(p,'t0',[])
addParameter(p,'domain',[]);

parse(p,obj,parmRange,varargin{:})
t0 = p.Results.t0;
objDomain = p.Results.domain;

%% Evaluate obj(s,t0) if specified
if isempty(t0)
    spaceParm = obj;
else
    t0 = varargin{1};
    spaceParm = obj.fixtime(t0);
end

spaceDimension = spaceParm.Dimension;
if size(parmRange,2) ~= 2*spaceDimension
    disp('we should stop here')
    error('parmRange should have 2*dimension-many columns')
end

%% Get spatial domain for parent surface
if nargout == 2
    if isempty(objDomain) % domain is [-1,1]^d
        spaceDomain = repmat([-1,1],1,spaceDimension);
    else
        spaceDomain = objDomain;
    end

    switch spaceDimension
        case 1
            % obtain new material coordinates
            rescale = @(s)(mean(spaceDomain) + .5*s*diff(spaceDomain));
            newMTC = rescale(parmRange);
        otherwise
            error('not yet implemented for higher dimensions')
    end
end

%% Reparameterize surface into subsurfaces
switch spaceDimension
    case 1 % spaceParm is a coefficient vector for a 1-d spatial parameterization.
        if length(spaceParm) ==1 % spaceParm is a single Scalar
            numSubSurface = size(parmRange,1); % number of subSurfaces
            if numSubSurface == 1 % parmRange specifies a single subSurface to be parameterized centered at (s1+s2)/2 with radius (s2-s1)/2.
                newCenter = intval(mean(parmRange)); % midpoint
                reScaleBy = intval(newCenter - parmRange(1)); % radius
                N = length(spaceParm.Coefficient);
                V = pascal(N);
                shiftOperator = intval(zeros(N));
                centerPowers = newCenter.^(0:N-1);
                for j = 1:N
                    shiftOperator(j,j:end) = V(j,1:N-j+1).*centerPowers(1:N-j+1)*reScaleBy^(j-1);
                end
                newParm = (shiftOperator*spaceParm.Coefficient')';


            else % each row of parmRange specifies a subsurface
                % if isequal(obj.NumericalClass,'intval')
                newParm = midrad(zeros(numSubSurface,spaceParm.Truncation),0);
                for j = 1:numSubSurface
                    newParm(j,:) = spaceParm.subdivide(parmRange(j,:));
                end
            end

        else % spaceParm is a vector of Scalars
            objLength = length(spaceParm);
            newParm{objLength} = spaceParm(objLength).subdivide(parmRange);
            for j = 1:objLength-1
                newParm{j} = spaceParm(j).subdivide(parmRange);
            end
        end

    otherwise % spaceParm is a higher dimension surface
        error('not yet implemented for higher dimensions')
end

switch nargout
    case 1
        varargout{1} = newParm; % no coordinate transform
    case 2
        varargout{1} = newMTC; % new material coordinates
        varargout{2} = newParm;
end
end % subdivide

% Revision History:
%{
14-Jun-2017 - implemented rigorous subdivision for Taylor parameterizations
13-Aug-2018 - updated for Scalar class
%}

