function varargout = eval(obj, data, varargin)
%EVAL - Evaluate a obj on space/time domain
%
%   Inputs:
%       obj - instance of Chart class
%       data - m-by-(k+1) of the form [S1,S2,...,Sk,T];
%
%   Outputs:
%       imageData - cell Array of evaluations in each coordinate
%       [Gamma_d,...,Gamma_d] - Evaluations for each scalar coordinate
%
%   Subfunctions: none
%   Classes required: @Chart, @Scalar
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jul-2018; Last revision: 18-Jul-2018

%%
% parse input and varargin
p = inputParser();
p.addRequired('obj')
p.addRequired('data')
p.addParameter('GlobalTime', false)
p.addParameter('GlobalSpace', false)
p.parse(obj, data, varargin{:})
globalTime = p.Results.GlobalTime;
globalSpace = p.Results.GlobalSpace;

% convert global to local coordinates if needed
if globalTime && globalSpace
    data = obj.intersectdomain(data);
    localTime = (data(:,end) - obj.TimeSpan(1))./obj.Tau; % convert global time to local time.
    localSpace = (2*data(:,1:end-1) - sum(obj.SpatialSpan))./diff(obj.SpatialSpan); % convert global time to local time.
    data = [localSpace, localTime];
    
elseif globalTime
    warning('this may omit endpoints due to floating point error')
    localTime = (data(:,end) - obj.TimeSpan(1))./obj.Tau; % convert global time to local time.
    validIdx = (0 <= localTime) & (localTime <= 1); % check data for evaluations which lie in this chart
    data = [data(validIdx, 1:end-1), localTime(validIdx)]; % filter out valid evaluations
    
elseif globalSpace
    warning('this may omit endpoints due to floating point error')
    
    localSpace = (2*data(:,1:end-1) - sum(obj.SpatialSpan))./diff(obj.SpatialSpan); % convert global time to local time.
    validIdx = (-1 <= localSpace) & (localSpace <= 1); % check data for evaluations which lie in this chart
    data = [localSpace(validIdx,:), data(validIdx,end)]; % filter out valid evaluations
end

%% evaluate data specified in local coordinates
if ~isempty(data)
    if nargout == obj.Dimension(2)
        for j = 1:obj.Dimension(2) % loop over phase space variables
            varargout{j} = real(obj.Coordinate(j).eval(data));
        end
    elseif nargout == obj.Dimension(2)+1
        for j = 1:obj.Dimension(2) % loop over phase space variables
            varargout{j} = real(obj.Coordinate(j).eval(data));
        end
        varargout{obj.Dimension(2)+1} = data;
    elseif nargout <= 1
        varargout{1} = {};
        for j = 1:obj.Dimension(2)
            varargout{1}{j} = real(obj.Coordinate(j).eval(data));
        end
    else
        error('nargout for eval all should be 1 or d')
    end
    
else
    if nargout == 1
        varargout{1} = {};
    else
        
        for j = 1:nargout
            varargout{j} = [];
        end
    end
end %  if
end %  eval

% Revision History:
%{
23-Mar-2019 - Added evaluation in global space/time coordinates
%}
