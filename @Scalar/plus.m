function sumObj = sum(obj)
%PLUS - Define addition for Scalar class
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

arrayCoefficient = [obj.Coefficient];
newCoefficient = reshape(arrayCoefficient,[],length(obj));
sumObj = Scalar(sum(newCoefficient,2)');
end % end plus

% Revision History:
%{
dd-Mon-yyyy - change1
dd-Mon-yyyy - change2
%}
