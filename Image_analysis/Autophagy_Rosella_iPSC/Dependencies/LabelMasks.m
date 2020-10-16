function [AutophagosomeMask, PhagophoreMask, LateAutolysosomeMask, EarlyAutolysosomeMask] = LabelMasks(ObjectList, Im)
%Creating category masks
%   ObjectList: Table containing class labels
%   Im: 3D image to use as a background

AutophagosomeObjects = ObjectList(ObjectList.AutophagosomeDecision == 1,:);
AutophagosomeMask = f_Create_Mask_from_ObjectList_Pixel_IDX(AutophagosomeObjects, 'PixelIdxList', Im);

PhagophoreObjects = ObjectList(ObjectList.PhagophoreDecision == 1,:);
PhagophoreMask = f_Create_Mask_from_ObjectList_Pixel_IDX(PhagophoreObjects, 'PixelIdxList', Im);

LateAutolysosomeObjects = ObjectList(ObjectList.LateAutolysosomeDecision == 1,:);
LateAutolysosomeMask = f_Create_Mask_from_ObjectList_Pixel_IDX(LateAutolysosomeObjects, 'PixelIdxList', Im);

EarlyAutolysosomeObjects = ObjectList(ObjectList.EarlyAutolysosomeDecision == 1,:);
EarlyAutolysosomeMask = f_Create_Mask_from_ObjectList_Pixel_IDX(EarlyAutolysosomeObjects, 'PixelIdxList', Im);


end

