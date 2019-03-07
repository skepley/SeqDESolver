classdef Chart < handle
    %CHART - A class for representing an analytic function, f: R^d ---> R^n.
    %
    %   CHART() - A more detailed description of the class
    %
    %   CHART constructor syntax:
    %       CHARTObj = $NAME()
    %
    %   CHART properties:
    %       Property 1 - description
    %       Property 2 - description
    %
    %   CHART methods:
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
    %   Other classes required: Scalar
    %   Other m-files required: none
    %   MAT-files required: none
    %
    %   Author: Shane Kepley
    %   email: shane.kepley@rutgers.edu
    %   Date: 21-Aug-2018; Last revision: 07-Mar-2019
    %
    %   ToDo:
    %   item1 -
    %   item2 -
    
    
    properties
        Coordinate; % length d vector of Scalars of dimension 1 + Dimension
        Dimension; % [d, n] 
        Truncation; % [N1,...,Nd]
        TimeSpan; % [t0, t1] - local coordinates for 1st dimension variable which is often time-like.
        SpatialSpan; % local coordinates with respect to [-1,1]^d for dimensions 2-d which are often space-like.
        Tau; % |t1 - t0|
        ErrorBound = 0; % rigorous error bound for validated charts or empty if not validated.
        FlowDirection = 1; % -1:backward time, 1:forward time
        NumericalClass; % double or intval
    end
    
    properties(Hidden = 1)
        Weight = 'ones'; % ell^1 space weights
        InitialData; % length n vector of Scalars of dimension d
        SubDivisionDepth; % recursion depth for use by subdivision algorithms
        InitialError = 0; % initial ell_1 error
        StepError = 0; % Validation error for this timestep (ErrorBound - InitialError)
        ParentHandle; % handle pointing to the parent chart if it exists
    end
    
end %end class

% Revision History:
%{
7-Mar-2019 - Re-factored for use with the new Atlas class.
%}

