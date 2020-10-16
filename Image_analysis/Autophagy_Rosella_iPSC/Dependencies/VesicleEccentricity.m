function [ObjectsPlusEccentricity] = VesicleEccentricity(ObjectList)
%Extract eccentricity by doing max projections cc by cc
%   Detailed explanation goes here


%% Extract Eccentricity by cc

EccentricityCell = cell(height(ObjectList), 1);

for c = 1:height(ObjectList)
    MaskThis3D = ObjectList.Image(c);
    MaskThis2D = max(MaskThis3D{:}, [], 3);
    CCThis = bwconncomp(MaskThis2D);
    ObjectThis = regionprops('table', CCThis, {'Eccentricity'});
    EccentricityCell{c,1} = ObjectThis.Eccentricity;
end

ObjectsPlusEccentricity = [ObjectList, EccentricityCell];
ObjectsPlusEccentricity.Properties.VariableNames{end} = 'Eccentricity';

end

