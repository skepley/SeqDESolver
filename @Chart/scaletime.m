function scaletime(chart, scaleBy)
%SCALETIME - Rescale time for a given chart based on single-step options

%   Syntax:
%       chart.SCALETIME(scaleBy)
%
%   Description:
%       SCALETIME() - 	% determine time scaling, L such that the time-1 map of the rescaled flow is equivalent to the time-L map of the given flow.
%
%   Inputs:
%       chart - instance of the Chart class
%       scaleBy - scalar double
%
%   Subfunctions: none
%   Classes required: @Chart, @Scalar
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Mar-2017; Last revision: 01-Jul-2018


%%
if ~strcmp(chart.Weight,'ones') % Assumption is that ell^1 weights should always be normalized to 1
    error('Weights other than one are not supported.')
end

if strcmp(chart.NumericalClass, 'intval')
    warning('Interval coefficients should never be rescaled this way.')
end
% 	if chart.FixTau == 0 % take a timestep based on coefficient decay but not to exceed MaxTau units
% 		rowDecay = chart.decay;
% 		scaleByDecay = (1e-16/rowDecay)^(1/(chart.Truncation(1)-1));
% 		scaleBy = min([scaleByDecay,chart.MaxTau]);
% 	else % take a fixed timestep
% 		scaleBy = chart.FixTau;
% 	end

chart.Coordinate.scaletime(scaleBy/chart.Tau) % call Scalar/scaletime on phase coordinates
chart.TimeSpan(2) = chart.TimeSpan(1) + chart.FlowDirection*scaleBy; % update chart timespan
chart.Tau = scaleBy; % update chart timestep

end %  scaletime
