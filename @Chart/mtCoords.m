%% EVALUATION METHODS
function MTC = mtCoords(obj,T)
    % An input real time vector returns a material time vector in [-1,1] relative to this timestep. 
	MTC = (T - obj.TimeSpan(1))/obj.Tau;
end