function varargout = local2global(obj, domainPt, varargin)
%CHART2REALTIME - Converts between chart (material) time and real time.
%
%   CHART2REALTIME() - A more detailed description of the function
%
%   Syntax:
%       tau = CHART2REALTIME(obj, t) evaluates the linear map:   tau(t):[-1,1] ---> [t0, tf] which maps chart time coordinates to actual times in the flow.
%       t = CHART2REALTIME(obj, tau, -1) evaluates the inverse linear map t(tau): [t0, tf] ---> [-1,1]
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

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Mar-2019; Last revision: 08-Mar-2019

if nargin == 3
    domainIdx = varargin{1};
else
    domainIdx = 1:length(domainPt);
end % end chart2realtime

if isequal(domainIdx, -1)
    error('Inverse mapping was removed from this method. This needs to be replaced by a global2local method.')
    %{
OLD CODE
tau = (2*domainPt(end) - sum(obj.TimeSpan))./obj.Tau; % linear interpolant for data [t0,tf], [-1,1] evaluated at t.
    %}
end

% map local (chart) coordinates to global (real) time coordinates
d = obj.Dimension(1); % get dimension of Chart
tau = []; % initialize as empty vectors
sigma = [];

if any(domainIdx < d) % some domain indices are spatial variables
    sigma = .5*(sum(obj.SpatialSpan) + domainPt(1)*diff(obj.SpatialSpan));
end

if ismember(d, domainIdx) % time variable is included in domainIdx
    tau = 0.5*(obj.TimeSpan(2)*(domainPt(end)+1) - obj.TimeSpan(1)*(domainPt(end)-1)); % linear interpolant for data [-1,1], [t0,tf] evaluated at t.
end
globalCoordinate = [sigma, tau]; % filter and return correct variable indices

switch nargout
    case 0 % just print to screen
        disp(globalCoordinate);
        
    case 1 % output is a d-vector of global [space, time] coordinates
        varargout{1} = globalCoordinate;
        
    case length(globalCoordinate) % unpack global coordinates
        varargout = mat2cell(globalCoordinate, 1, length(globalCoordinate));
        
    otherwise
        error('Number of output arguments must be 0, 1, or d')
end


% Revision History:
%{
08-Mar-2019 - This file replaces the previous version mtCoords.
17-Jun-2019 - Updated documentation. Fixed output for single varargout calls. Removed inverse mapping function. Added support for specifying
    partial indices to return global domain values.
%}
