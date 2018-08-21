classdef Chart < handle
%CHART - A class for representing an analytic function, f: R^d ---> R^n. 
    
    properties
		Coordinate; % length d vector of Scalars of dimension 1 + Dimension
        Dimension; % [d, n]
        Truncation; % [N1,...,Nd]
        TimeSpan; % [t0, t1]
        SpatialSpan; % material coordinates with respect to [1,1]^d
		Tau; % |t1 - t0|
        ErrorBound = 0; % double or empty array
        FlowDirection = 1; % -1:backward time, 1:forward time
        NumericalClass; % double or intval
    end

    properties(Hidden = 1)
		Weight = 'ones'; % ell^1 space weights
        InitialData; % length n vector of Scalars of dimension d
        SubDivisionDepth; % recursion depth for use by subdivision algorithms
        InitialError = 0; % initial ell_1 error
        StepError = 0; % Validation error for this timestep (ErrorBound - InitialError)
        ParentHandle; % handle pointing to the parent chart
    end

end %end class


