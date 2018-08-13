function varargout = eval(chart,data,varargin)
%EVAL - Evaluate a chart on space/time domain 
%
%   Inputs:
%       chart - instance of Chart class
%       data - m-by-(k+1) of the form [S1,S2,...,Sk,T];
%
%   Outputs:
%       imageData - cell Array of evaluations in each coordinate
%       [Gamma_d,...,Gamma_d] - Evaluations for each scalar coordinate
%
%   Subfunctions: none
%   Classes required: @Chart, @BAscalar
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jul-2018; Last revision: 18-Jul-2018

%%
% parse input and varargin
p = inputParser();
p.addRequired('chart')
p.addRequired('data')
p.addParameter('AbsoluteTime',false)
p.addParameter('AbsoluteSpace',false)

p.parse(chart,data,varargin{:})
absoluteTime = p.Results.AbsoluteTime;
absoluteSpace = p.Results.AbsoluteSpace;


%% evaluate data specified in local coordinates
if ~isempty(data)
    if nargout == chart.PhaseDimension
        for j =1:chart.PhaseDimension
            varargout{j} = real(chart.Coordinate(j).eval(data));
        end
        
    elseif nargout ==1
        varargout{1} = {};
        for j =1:chart.PhaseDimension
            varargout{1}{j} = real(chart.Coordinate(j).eval(data));
        end
    else
        error('nargout for eval all should be 1 or d')
    end
    
else
    for j = 1:nargout
        varargout{j} = [];
    end
end % end if
end % end eval

