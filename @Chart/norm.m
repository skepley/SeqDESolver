function BAnorm = norm(obj)
    phaseNorms = obj.var.norm;
    DFNorms = obj.DFvar.norm;
    otherNorms = obj.othervar.norm;
    BAnorm = max([max(phaseNorms(:)),max(DFNorms(:)),max(otherNorms(:))]);
    obj.Norm = BAnorm;
end