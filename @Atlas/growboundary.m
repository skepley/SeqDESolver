function [output1,output2] = growboundary(boundaryChart, subdividecheck, advectioncheck, evaluationcheck, varargin)
%GROWBOUNDARY - A method for growing Atlas by integrating the boundary
%
%   GROWBOUNDARY() - Performs a single iteration of the following 3 steps: subdivision, advection, and evaluation. This takes k^th generation boundary 
%   leaves in the chart tree to the (k+1)^st generation boundary by advecting.
%
%   Syntax:
%       output = GROWBOUNDARY(input1, input2)
%       [output1, output2] = GROWBOUNDARY(input1, input2, input3)
%    
%   Inputs:
%       boundaryCharts - Description
%       boundarycheck - A function which decides if boundary charts should be passed to subdivision phase. If not they crash.
%       subdividewhere - A function which returns spatial nodes where a boundary chart should be subdivided.
%       advectioncheck - A function which decides if boundary charts should be advected. If not they crash.
%       evaluationcheck - A function which decides if advected charts should be evaluated to boundaries. If not they crash.

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
%   Date: 07-Mar-2019; Last revision: 07-Mar-2019

%% parse input 
p = inputParser;
addRequired(p, boundaryChart)
addRequired(p, boundarycheck)
addRequired(p, subdividewhere)
addRequired(p, advectioncheck)
addRequired(p, evaluationcheck)

addParameter(p,'Parameter1',default1)
addParameter(p,'Parameter2',default2)
addParameter(p,'Parameter3',default3)

parse(p,boundaryChart,subdividecheck,input3,varargin{:})
parameter1 = p.Results.Parameter1;
parameter2 = p.Results.Parameter2;
parameter3 = p.Results.Parameter3;


%% ---------------------------------------- SUBDIVISION PHASE ----------------------------------------
for j = 1:length(boundaryChart)
   if boundarycheck(boundaryChart(j)) && subdividecheck(boundaryChart(j))
       jSubDivisionNode = subdividewhere(boundaryChart(j)); % get nodes for subdivision
   else
       
   end
    
end


%% ---------------------------------------- ADVECTION PHASE ----------------------------------------


%% ---------------------------------------- EVALUATION PHASE ----------------------------------------



end % end growboundary

% Revision History:
%{

%}
