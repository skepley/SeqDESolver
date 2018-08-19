function varargout = ellonebox(chart)
%ELLONEBOX - returns a 1-by-n interval vector which bounds the chart in phase space using coarse ell^1 bounds.

%   Syntax:
%       chart.ELLONEBOX() - Append ell^1 box to chart.ellOneBox property
%       box = chart.ELLONEBOX() - Append and return box
%
%   Description:
%       ELLONEBOX() - description
%
%   Inputs:
%       chart - instance of the Chart class
%
%   Outputs:
%       box - 1-by-n intval
%
%   Other m-files required: @Chart, @Scalar
%   Subfunctions: Scalar/norm
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 01-Jul-2018; Last revision: 01-Jul-2018

if length(chart) > 1 % vectorized version
    for j = 1:length(chart)
        chart(j).ellonebox();
    end

elseif ~isempty(chart.ellOneBox) % ellOneBox is already stored
    box = chart.ellOneBox;

else % computer box for a single chart
    center = arrayfun(@(j)chart.Coord(j).Coef(1,1), 1:chart.Dimension(2)); % constant term in chart expansion
    fullNorm = chart.Coord.norm;
    radius = fullNorm - abs(center) + chart.StepError; % chart.StepError = 0 for nonrigorous timesteps.
    box = midrad(center, radius)';
    chart.ellOneBox = box;
end

varargout{1} = box; % return box if argout
end

