function [sigma, tau] = local2global(obj, domainPt, varargin)
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
    transformDirection = varargin{1};
else
    transformDirection = 1; % default
end % end chart2realtime

switch transformDirection
    case 1 % map local (chart) coordinates to global (real) time coordinates
        sigma = .5*(sum(obj.SpatialSpan) + domainPt(1)*diff(obj.SpatialSpan));
        tau = 0.5*(obj.TimeSpan(2)*(domainPt(end)+1) - obj.TimeSpan(1)*(domainPt(end)-1)); % linear interpolant for data [-1,1], [t0,tf] evaluated at t.
    case -1 % compute the inverse map taking real time to chart time
        tau = (2*domainPt(end) - sum(obj.TimeSpan))./obj.Tau; % linear interpolant for data [t0,tf], [-1,1] evaluated at t.
        error('need formula for spatial map from global to local')
end
% Revision History:
%{
08-Mar-2019 - This file replaces the previous version mtCoords.
%}
