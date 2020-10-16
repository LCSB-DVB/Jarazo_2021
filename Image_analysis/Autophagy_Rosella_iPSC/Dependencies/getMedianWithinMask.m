function [TheMedian] = getMedianWithinMask(ImInCell, MaskinCell)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    Im = ImInCell{:};
    Mask = logical(MaskinCell{:});
    TheMedian = median(Im(Mask));

end

