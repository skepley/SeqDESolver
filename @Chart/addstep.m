
        
        %% TIMESCALING AND TIMESTEPPING METHODS
        function tstep = step(obj,tau,varargin)
            coefs_t = obj.tau_map(tau);
            switch nargin
                case 2
                    tstep = vdp_timestep(obj.params,coefs_t(1,:),coefs_t(2,:),obj.modes,'tspan',obj.tspan(2));
                case 3
                    forcing = varargin{1};
                    tstep = vdp_timestep(obj.params,coefs_t(1,:),coefs_t(2,:),obj.modes,'tspan',obj.tspan(2),'forcing',forcing);
                case 4
                    forcing = varargin{1};
                    L = varargin{2};
                    tstep = vdp_timestep(obj.params,coefs_t(1,:),coefs_t(2,:),obj.modes,'tspan',[obj.tspan(2), obj.tspan(2) + L],'forcing',forcing);
            end
        end