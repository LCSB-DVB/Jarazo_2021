function [PerimMask3D] = f_bwperim3D(Mask3D)
%Create perimeter mask in 3D
%Author: Paul Antony 2016/11/04 
%   Mask3D: 3 dimensional input mask

PerimMask3D = Mask3D - imerode(Mask3D, strel('sphere',1));

end

