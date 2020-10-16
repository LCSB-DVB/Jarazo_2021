function [LabelMask] = LabelCategories(ObjectList, Im)
%Creating a Mask with all the categories in differnt colors shown on each
%plane of 3D Im
%   ObjectList: Table containing class labels
%   Im: 3D image to use as a background

AutophagosomeObjects = ObjectList(ObjectList.AutophagosomeDecision == 1,:);
AutophagosomeMask = f_Create_Mask_from_ObjectList_Pixel_IDX(AutophagosomeObjects, 'PixelIdxList', Im);
AutophagosomeMaskPerim = f_bwperim3D(AutophagosomeMask);

PhagophoreObjects = ObjectList(ObjectList.PhagophoreDecision == 1,:);
PhagophoreMask = f_Create_Mask_from_ObjectList_Pixel_IDX(PhagophoreObjects, 'PixelIdxList', Im);
PhagophoreMaskPerim = f_bwperim3D(PhagophoreMask);
PhagophoreMaskPerim = PhagophoreMaskPerim .* 2;

LateAutolysosomeObjects = ObjectList(ObjectList.LateAutolysosomeDecision == 1,:);
LateAutolysosomeMask = f_Create_Mask_from_ObjectList_Pixel_IDX(LateAutolysosomeObjects, 'PixelIdxList', Im);
LateAutolysosomeMaskPerim = f_bwperim3D(LateAutolysosomeMask);
LateAutolysosomeMaskPerim = LateAutolysosomeMaskPerim .* 3;

EarlyAutolysosomeObjects = ObjectList(ObjectList.EarlyAutolysosomeDecision == 1,:);
EarlyAutolysosomeMask = f_Create_Mask_from_ObjectList_Pixel_IDX(EarlyAutolysosomeObjects, 'PixelIdxList', Im);
EarlyAutolysosomeMaskPerim = f_bwperim3D(EarlyAutolysosomeMask);
EarlyAutolysosomeMaskPerim = EarlyAutolysosomeMaskPerim .* 4;

ClassifiedMask = AutophagosomeMaskPerim + PhagophoreMaskPerim + LateAutolysosomeMaskPerim + EarlyAutolysosomeMaskPerim;

LabelMask = f_imoverlay3D_5Masks(Im, AutophagosomeMaskPerim, PhagophoreMaskPerim, LateAutolysosomeMaskPerim, EarlyAutolysosomeMaskPerim, zeros(size(Im), 'logical'), {'AutophagosomeMaskPerim\_R PhagophoreMaskPerim\_G LateAutolysosomeMaskPerim\_B EarlyAutolysosomeMaskPerim\_V'})
% add saving to the function f_imoverlay3D_5Masks
end

