function setinitialdata(chart,initialData)
%SETINITIALDATA - Write initial data to chart as a BAscalar
%
%
%   Inputs:
%       chart - instance of Chart class
%       initialData - (1+d)-dimensional array. (j,:,...,:) is the j^th coordinate in phase space defined as a d-dimsional coefficient array.
%
%   Outputs:
%       chart - chart instance with initial data appended
%
%   Other m-files required: @Chart, @BAscalar
%   Subfunctions: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 01-Jul-2018; Last revision: 01-Jul-2018


% append initial data to chart
chart.InitialData = BAscalar();
switch chart.SurfaceDimension
    case 0 % initial data is a point
        for j = 1:chart.PhaseDimension
            chart.InitialData(j) = BAscalar(initialData(j,:));
        end
        
    case 1 % initial data is an arc
        for j = 1:chart.PhaseDimension
            chart.InitialData(j) = BAscalar(initialData(j,:), chart.SpatialTruncation);
        end
    case 2 % initial data is a 2-d surface
        for j = 1:chart.PhaseDimension
            chart.InitialData(j) = BAscalar(initialData(j,:,:), chart.SpatialTruncation);
        end
    otherwise
        error('not implemented for higher dimensions')
end
end

