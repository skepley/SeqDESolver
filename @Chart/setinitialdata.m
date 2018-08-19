function setinitialdata(obj, initialData, basis)
%SETINITIALDATA - Write initial data to chart as a Scalar
%
%   Inputs:
%       chart - instance of Chart class
%       initialData - (1+d)-dimensional array. (j,:,...,:) is the j^th coordinate in phase space defined as a d-dimsional coefficient array.
%
%   Outputs:
%       chart - chart instance with initial data appended
%
%   Other m-files required: @Chart, @Scalar
%   Subfunctions: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 01-Jul-2018; Last revision: 14-Aug-2018

% append initial data to chart
obj.InitialData = Scalar();
switch obj.Dimension(1)-1
    case 0 % initial data is a point
        for j = 1:obj.Dimension(2)
            obj.InitialData(j) = Scalar(initialData(j,:), basis);
        end
    case 1 % initial data is an arc
        for j = 1:obj.Dimension(2)
            obj.InitialData(j) = Scalar(initialData(j,:), basis, obj.Truncation(2:end));
        end
    case 2 % initial data is a 2-d surface
        for j = 1:obj.Dimension(2)
            obj.InitialData(j) = Scalar(initialData(j,:,:), basis, obj.SpatialTruncation(2:end));
        end
    otherwise
        error('not implemented for higher dimensions')
end
end

% Revision History:
%{
14-Aug-2018 - updated for Scalar class
%}
