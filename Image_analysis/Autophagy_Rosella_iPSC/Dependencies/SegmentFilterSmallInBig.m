function [MaskToAdd] = SegmentFilterSmallInBig(MaskBig, MaskSmall, SizeThreshold, ProportionThreshold)
%This function merges small connected components with bigger touching
%connected components in 2D or 3D
%   MaskBig: the mask to be processed
%   MaskSmall: the mask of contained vesicles
%   SizeThreshold: Threshold separating small and big connected components
%   ProportionThreshold: Volume proportion
%   MaskToAdd: The big objects or pixels within big objects to be kept as seeds
    
    MaskToAdd = zeros(size(MaskBig));
    BigMask = bwareaopen(MaskBig, SizeThreshold); %vol(BigMask)

    BigCC = bwconncomp(BigMask);
    BigObjects = regionprops('table', BigCC, {'PixelIdxList'});
    
    if height(BigObjects) == 0
        MaskReconstructed = zeros(size(MaskBig), 'logical');
        return
    end
    
    for o = 1:height(BigObjects)
        BigObjectThis = BigObjects(o,:);
        MaskThis = zeros(size(MaskBig), 'logical');
        MaskThis = f_Create_Mask_from_ObjectList_Pixel_IDX(BigObjectThis, 'PixelIdxList', MaskBig); % vi(uint8(MaskThis))
        
        %%
        ContainedVesicleMask = MaskThis .* MaskSmall; % vol(ContainedVesicleMask)
        ContainedVesicleMetric = sum(ContainedVesicleMask(:)) / sum(MaskThis(:));
        
        if ContainedVesicleMetric > ProportionThreshold
            MaskToAdd = MaskToAdd + MaskThis .* MaskSmall;
        else
            MaskToAdd = MaskToAdd + MaskThis;
        end
        
        
    end

    MaskToAdd  = logical(MaskToAdd );
    %vol(MaskToAdd )

end

