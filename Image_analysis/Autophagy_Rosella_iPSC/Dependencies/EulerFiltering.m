function [ MaskOut ] = EulerFiltering( EulerIm )
%Mask and filter a contour image coming from fourrier filtering to identify
%objects with one autophagosome alike hole
%   Detailed explanation goes here

    EulerMask = EulerIm > 150;
    %vi(double(flip(EulerIm,1)))
    %imwrite(imadjust(uint16(EulerMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\GreenDoG150.png')
    
    %vi(double(flip(EulerMask,1)))
    
    EulerMaskRings3D = EulerMask; % To use later to extract relevant planes
    EulerMask = max(EulerMask,[],3);
    EulerMask = bwareaopen(EulerMask, 20);
    EulerMask = medfilt2(EulerMask);
    EulerMask = imopen(EulerMask, strel('disk',1));
    %imwrite(imadjust(uint16(EulerMask)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\GreenDoG150b.png')
    
    % filter based on Euler number
    EulerCC = bwconncomp(EulerMask);
    EulerObjects = regionprops('table', EulerCC, {'EulerNumber','Eccentricity','PixelIdxList'});
    EulerObjects = EulerObjects(EulerObjects.EulerNumber == 0,:);
    EulerMask = f_Create_Mask_from_ObjectList_Pixel_IDX(EulerObjects, 'PixelIdxList', EulerMask);
    %imwrite(imadjust(uint16(EulerMask)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\EulerZero.png')
    
    % filter good candidates with expected concavity proportion and size
    EulerMaskFilled = imfill(EulerMask);
    %imwrite(imadjust(uint16(EulerMaskFilled)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\EulerZeroFilled.png')

    EulerMask_CC = bwconncomp(EulerMask);
    EulerMaskFilled_Objects = regionprops('table', EulerMask_CC,{'PixelIdxList', 'Area', 'FilledArea'});
    EulerMaskFilled_Objects = EulerMaskFilled_Objects((EulerMaskFilled_Objects.FilledArea ./ EulerMaskFilled_Objects.Area) > 1.01, :);
    EulerMaskFilled_Objects = EulerMaskFilled_Objects((EulerMaskFilled_Objects.FilledArea - EulerMaskFilled_Objects.Area) > 20, :);
    EulerMask = f_Create_Mask_from_ObjectList_Pixel_IDX(EulerMaskFilled_Objects, 'PixelIdxList', EulerMask);
    %imwrite(imadjust(uint16(EulerMask)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\EulerSelect.png')
    
    % Create the actual vesicle mask
    EulerMask = ~EulerMask; % vol(EulerMask,0,1)
    EulerMask = f_RemoveBigObjects(EulerMask, 1e4);
    EulerMask = imreconstruct(EulerMask, EulerMaskFilled);
    EulerMask = imopen(EulerMask, strel('disk',5));
    EulerMask = cat(3, EulerMask, EulerMask, EulerMask, EulerMask, EulerMask);
    %imwrite(imadjust(uint16(EulerMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\EulerVesicles.png')

    %% Keep most relevant planes of each connected component
    
    EulerMask_CC2 = bwconncomp(EulerMask);
    EulerMask2_Objects = regionprops('table', EulerMask_CC2, EulerMaskRings3D,{'PixelIdxList', 'Image', 'Area', 'SubarrayIdx'}); % Every objects still uses all planes
    PlaneWeightMatrix = [];

    for o = 1:height(EulerMask2_Objects) %% loop over connected components
        CCim = EulerMask2_Objects(o,'Image'); CCim = CCim{1,1}{:};

        Idx = EulerMask2_Objects(o,'SubarrayIdx'); Idx = Idx{1,:};
        WeightCCim = EulerMaskRings3D(Idx{1},Idx{2},Idx{3});

        for p = 1:5 % Loop over planes and collect plane by plane features in a connected component table
            CCimPlane = CCim(:,:,p);
            WeightCCimPlane = WeightCCim(:,:,p);
            PlaneWeightMatrix(o,p,1) = sum(CCimPlane(:)); % Plane filled
            PlaneWeightMatrix(o,p,2) = sum(WeightCCimPlane(:)); % Plane from RingMask
        end

        PlaneWeightMatrix(:,:,3) = PlaneWeightMatrix(:,:,2) ./ PlaneWeightMatrix(:,:,1); % Normalized metric
        PlaneWeightMatrix(:,:,4) = PlaneWeightMatrix(:,:,3) > 0.25; % Decision

    end
    
    if height(EulerMask2_Objects) == 0 % handle the situation of images without any retained vesicles for the output mask
        
        MaskOut = zeros(size(EulerIm), 'logical');
        
    else

        for o = 1:height(EulerMask2_Objects) % add retained vesicles to the output mask

            ObjectThis = EulerMask2_Objects(o,:);
            MaskThis = f_Create_Mask_from_ObjectList_Pixel_IDX(ObjectThis, 'PixelIdxList', EulerMask);
            MaskThis(:,:,~logical(PlaneWeightMatrix(o,:,4))) = 0; % removing planes with low ring mask contribution

            if o==1
                MaskOut = MaskThis;
            else
                MaskOut = MaskThis | MaskOut;
            end

        end

        MaskOut = imopen(MaskOut, strel('sphere',1)); % smoothing in 3D
    
    end
    

end

