classdef Atlas
    %ATLAS - A class for working with collections of Charts.
    %
    %   ATLAS() - A more detailed description of the class
    %
    %   ATLAS constructor syntax:
    %       AtlasObj = $NAME()
    %
    %   ATLAS properties:
    %       Property 1 - description
    %       Property 2 - description
    %
    %   ATLAS methods:
    %       Method 1 - description
    %       Method 2 - description
    %
    %   Examples:
    %       Line 1 of example
    %       Line 2 of example
    %       Line 3 of example
    %
    %   Subclasses: none
    %   Superclasses: none
    %   Other classes required: Chart, Scalar
    %   Other m-files required: none
    %   MAT-files required: none
    %
    %   Author: Shane Kepley
    %   email: shane.kepley@rutgers.edu
    %   Date: 07-Mar-2019; Last revision: 07-Mar-2019
    %
    %   ToDo:
    %   item1 -
    %   item2 -
    
    %% -------------------- Properties --------------------
    properties
        ChartTree; % a tree of Chart objects ordered linearly
        Size; % number of Charts in the Atlas
        Leaf; % index for leaves in the tree
        Height; % height of the ChartTree which is equivalent to the maximum number of generations
        Valid; % logical for whether the charts are rigorously validated or not
        FlowDirection = 1; % -1:backward time, 1:forward time
    end
    
    properties(Hidden = 1)
        SubDivisionOptions; % options for subdivision algorithms
    end
    
    
    %% -------------------- Methods --------------------
    methods
        function obj = Atlas(boundaryChart,input2,input3,varargin)
            %ATLAS - class constructor
            %
            % initialization should be a single boundary chart
            
        end % end class constructor
    end % end methods
    
end % end classdef



% Revision History:
%{

%}