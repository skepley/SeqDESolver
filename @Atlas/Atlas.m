classdef Atlas < handle
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
        Chart; % a tree of Chart objects ordered linearly
        Size=0; % number of Charts in the Atlas
        LeafStack; % a vector of handles to the boundary charts which are leaves of the tree
        CrashStack; % a vector of handles for boundary charts which have crashed (stopped integration) for any reason.
        Height; % height of the ChartTree which is equivalent to the maximum number of generations
        Valid; % logical for whether the charts are rigorously validated or not
%         FlowDirection = 1; % -1:backward time, 1:forward time
        TimeStepper; % handle for function which decides how to timestep
        Truncator; % handle for function which decides how to truncate
        BoundaryCheck % handle for checking whether to pass a boundary Chart to advection stage
        AdvectionCheck % handle for checking whether to accept an advected Chart
        MaxTau; % temporary property
        MaxGeneration;
    end
    
    properties(Hidden = 1)
        SubDivisionOptions; % options for subdivision algorithms
        ChartClass; % name of Chart class for this atlas
        LastGeneration = 0;
    end
    
    
    %% -------------------- Methods --------------------
    methods
        function obj = Atlas(boundaryChart, tau, boundarycheck, advectioncheck, varargin)
            %ATLAS - class constructor

            if nargin > 0
                % initialization should be a single boundary chart
                
                % parse input
                p = inputParser;
                addRequired(p,'boundaryChart') 
                addRequired(p,'tau') 
                addRequired(p,'boundarycheck') % 1-by-1 double
                addRequired(p,'advectioncheck') % 1-by-(1+d) integer
                addParameter(p,'MaxTau', Inf) % scalar double
                addParameter(p,'MaxGeneration', Inf) % logical
                
                % parse variable arguments
                parse(p, boundaryChart, tau, boundarycheck, advectioncheck, varargin{:})
                boundaryChart = p.Results.boundaryChart; % set the regularization type
                tau = p.Results.tau;
                obj.MaxTau = p.Results.MaxTau;
                obj.MaxGeneration = p.Results.MaxGeneration; % generate the initial data as a boundary chart but don't advect it.
                obj.LeafStack = boundaryChart;
                obj.TimeStepper = @(bdChart)tau; % fixed timestep
                obj.ChartClass = boundaryChart(1).ClassConstructor();
                obj.CrashStack = obj.ChartClass; % initialize crash stack
                obj.CrashStack = obj.CrashStack(2:end); % empty the stack
                obj.BoundaryCheck = p.Results.boundarycheck; 
                obj.AdvectionCheck = p.Results.advectioncheck;
            end % end if nargin
        end % end class constructor
    end % end methods
    
end % end classdef



% Revision History:
%{

%}