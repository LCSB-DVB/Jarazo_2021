function [Mask3Dfiltered] = f_RemoveNonSphericalObjects(Mask3D, SphericityThreshold1, SphericityThreshold2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Mask3DPadder = zeros([size(Mask3D, 1), size(Mask3D, 2), size(Mask3D, 3)+2], 'uint16');
Mask3DPadder(:,:,2:end-1) = Mask3D;
Mask3DPadder = imclearborder(Mask3DPadder); % Objects touching top and bottom are protected from added zero planes
Mask3D = Mask3DPadder(:,:,2:end-1);

cc = bwconncomp(Mask3D);
% Objects = regionprops('table', cc, {'Area', 'FilledArea','PixelIdxList','MaxIntensity'});
% Objects.Properties.VariableNames{4} = 'MinimumFound';

% Surface analysis          
ErodedMask = imerode(Mask3D, strel('sphere',1));
SurfMask = (Mask3D - ErodedMask) > 0;
%vi(flip(uint8(SurfMask),1))

Objects = regionprops('table', cc, SurfMask,{'Area', 'PixelValues', 'PixelIdxList'});
SurfPixelCounts = rowfun(@(x) sum(x{:}), Objects, 'InputVariables', 'PixelValues');
%vi(uint8(SurfMask))

Objects = [Objects, SurfPixelCounts];
Objects.Properties.VariableNames{end} = 'Surface';



SphericityIndex1 = rowfun(@(x,y) ((pi^(1/3))*((6*x)^(2/3)))/y, Objects, 'InputVariables', {'Area', 'Surface'}); % Higher values indicate better sphericity
SphericityIndex1.Properties.VariableNames = {'SphericityIdx1'};
SphericityIndex2 = rowfun(@(x,y) x/y, Objects, 'InputVariables', {'Area', 'Surface'}); % Higher values indicate better sphericity
SphericityIndex2.Properties.VariableNames = {'SphericityIdx2'};

Objects = [Objects, SphericityIndex1, SphericityIndex2];
%Objects.Properties.VariableNames{end-1:end} = {'SphericityIdx1', 'SphericityIdx2'};

ObjectsSpherical = Objects(Objects.SphericityIdx1 > SphericityThreshold1 & Objects.SphericityIdx2 > SphericityThreshold2,:);
Mask3Dfiltered = f_Create_Mask_from_ObjectList_Pixel_IDX(ObjectsSpherical, 'PixelIdxList', Mask3D);

% disp('debug')
end

