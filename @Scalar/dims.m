function objDims = dims(obj)
% returns the dimensions for the Scalar
objDims = size(obj.Coefficient) > 1;
end
