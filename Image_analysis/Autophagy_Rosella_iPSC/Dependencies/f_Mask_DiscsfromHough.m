function [MaskDisks, MaskSplitDisks] = f_Mask_DiscsfromHough(Im, HTG_centers, HTG_radii, Candidates)
%Creates a Mask of discs found via a HoughTransform (see imfindcircles)
%   Paul Antony 20161103
%   Im: image of desired output size (as used previously for hough function)
%   HTG_centers: xy coordinates of disc centers from Hough
%   HTG_radii: radii
%   Candidates: boolean or index vector of discs to add to mask

CircleBackground = zeros(size(Im), 'uint8');
CirclesFound = [HTG_centers, HTG_radii];
CirclesToKeep = CirclesFound(Candidates,:);
MaskDisks = insertShape(CircleBackground, 'FilledCircle', CirclesToKeep, 'SmoothEdges', false);
MaskDisks = max(MaskDisks,[],3) > 0; %it(MaskDisks)

MaskDisksCenters = CircleBackground;
Idx = round(HTG_centers);
IdxLin = sub2ind(size(Im), Idx(Candidates,2), Idx(Candidates,1));
MaskDisksCenters(IdxLin) = 1; % it(MaskDisksCenters)

WS_MaskDisksCenters = bwdist(MaskDisksCenters);%it(WS_MaskDisksCenters)
MaskSplitDisks = watershed(WS_MaskDisksCenters);
MaskSplitDisks = uint16(MaskSplitDisks) .* uint16(MaskDisks);
MaskSplitDisks = imopen(MaskSplitDisks, strel('disk', 2));%it(MaskSplitDisks)
MaskSplitDisks = logical(MaskSplitDisks);


end

