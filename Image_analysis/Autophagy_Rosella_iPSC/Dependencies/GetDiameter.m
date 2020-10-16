function [Diameter] = GetDiameter(Mask)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Mask = Mask{:};
Mask = max(Mask,[],3);
CC = bwconncomp(Mask);
Object = regionprops('table', CC, {'MajorAxisLength'});
Diameter = table2array(Object);

end

