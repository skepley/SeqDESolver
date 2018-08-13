classdef Chart < handle
    % A vector of k-dimensional BAscalars of length d. timestep(i) represents a timestep of the i^th coordinate in R^d for
    
    properties
		Coordinate; % length d vector of BAscalars of dimension 1 + SurfaceDimension
        SpatialDegree; % degree of spatial polynomial truncation
        TemporalDegree; % degree of temporal polynomial truncation
        TimeSpan; % [t0, t1]
        SpatialSpan; % material coordinates with respect to [1,1]^d
		Tau; % |t1 - t0|
        ErrorBound = 0; % double or empty array
        FlowDirection = 1; % -1:backward time, 1:forward time
        CoefType; % double or intval
        StepError = 0; % Validation error for this timestep
    end
    
    properties(Hidden = 1)
        PhaseDimension; % Dimension of phase space
		SurfaceDimension; % Dimension of initial data
		Weight = 'ones'; % ell^1 space weights
        InitialData; % length d vector of BAscalars of dimension SurfaceDimension
		ErrorProp;
        SubArcDepth;
        InitialError = 0;
        ParentHandle; % handle pointing to the parent chart
    end
    
    % methods
        % function obj = taylor_timestep(coef_gen,params,init_surface,modes,varargin)
            % % class constructor
			% % specify init_surface as a d-vector of k-dimensional BAscalars
            % if(nargin > 0)
                % p = inputParser;
				% addRequired(p,'coef_gen')
                % addRequired(p,'params')
                % addRequired(p,'init_surface')
                % addRequired(p,'modes')
                % addRequired(p,'tspan'); % [t0,t1] for fixed timestep, or just t0 for longest timestep with autorescaling
                % addParameter(p,'forcing',{0,0});

                % % parse varargs
                % parse(p,params,init_surface,modes,varargin{:})
                % obj.forcing = p.forcing;
                
                % % set properties
                % obj.params = params;
                % obj.modes = modes;
				% obj.phase_dimension = length(init_surface);
				% obj.surface_dimension = init_surface(1).surface_dimension;
				% coefs_gen(init_surface);
                % switch obj.surface_dimension
                    % case 0 % initial data is a point
                        % % obj.var = [BAscalar(X_init);BAscalar(Y_init)];
						% % No scaling on IVP solver. 
                        
                    % case 1 % initial data is an arc
                        % % obj.var = [BAscalar(X_init,[1,obj.modes(2)]);BAscalar(Y_init,[1,obj.modes(2)])];
                        % if length(obj.tspan) == 1
                            % auto_scale(obj);
                        % else
                            % obj.scale_time(diff(obj.tspan));
                        % end
                % end
            % end
        % end
    % end %end methods
end %end class


