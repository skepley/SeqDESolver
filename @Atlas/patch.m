function varargout = patch(obj, evalNode, idx, varargin)
%PATCH - Plot an atlas of Charts as a patch
%
%   PATCH() - A more detailed description of the function
%
%   Syntax:
%       output = PATCH(input1, input2)
%       [output1, output2] = PATCH(input1, input2, input3)
%    
%   Inputs:
%       obj - An atlas of Chart objects parameterized on [-1,1] x [0, Tau]
%       evalNodes - A cell array of global space and time evaluation nodes
%       idx - coordinate indices to plot
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
%   Date: 19-Apr-2019; Last revision: 19-Apr-2019

%% parse input 
% p = inputParser;
% addRequired(p,evalNode)
% addRequired(p,input2)
% addRequired(p,input3)
% addParameter(p,'Parameter1',default1)
% addParameter(p,'Parameter2',default2)
% addParameter(p,'Parameter3',default3)
% 
% parse(p,evalNode,input2,input3,varargin{:})
% parameter1 = p.Results.Parameter1;
% parameter2 = p.Results.Parameter2;
% parameter3 = p.Results.Parameter3;

atlasDimension = length(evalNode);

% spaceNode = linspace(-1,1,500);
% timeNode = linspace(0,.33,500);
switch atlasDimension
    case 2
        spaceNode = evalNode{1};
        timeNode = evalNode{2};
        [S,T] = meshgrid(spaceNode,timeNode);
        evalData = [reshape(S,[],1), reshape(T,[],1)];
    otherwise
        error('not implemented')
end

% initialize patch data
face = [];
vertex = [];
colorData = [];
% main loop
for iChart = obj.Chart
    % get global evaluation node data for this chart
    bdVertices = [iChart.SpatialSpan(1), iChart.TimeSpan(1); iChart.SpatialSpan(1), iChart.TimeSpan(2);... % add 4 corners of domain
        iChart.SpatialSpan(2), iChart.TimeSpan(1); iChart.SpatialSpan(2), iChart.TimeSpan(2)];
    data = [iChart.intersectdomain(evalData);bdVertices];
    
    % color data is time
    % iColorData = data(:,2);
    colorData = [colorData; data(:,2)];
    

    % get global vertex data for this chart
    iPatchCell = iChart.eval(data, 'GlobalTime', true, 'GlobalSpace', true); % cell array of evaluations
    iPatchEval = cell2mat(iPatchCell(idx)); % convert to data array consisting only of indices to plot
    
    % jVertex = [mid(u),mid(v),mid(q)];
    iVertex = mid(iPatchEval); % convert to floats if necessary
    vertexCount = size(vertex, 1); % current vertex count to shift global vertex indices of this chart
    vertex = [vertex;iVertex]; % append new vertices
    
    % get global face data for this chart
    switch atlasDimension
        case 2
            iFace = delaunay(data(:,1),data(:,2)); % index triples for delaunay triangulation in space-time.
        otherwise
            error('not implemented')
    end
    face = [face;iFace + vertexCount]; % start counting face labels from the previous end label
end

% plot it
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none');
end % end patch

% Revision History:
%{

%}
