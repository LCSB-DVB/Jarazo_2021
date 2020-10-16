function [ RawIm3D ] = ApplyFlatFieldCorrection(RawIm3D, CorrIm2D)
%Apply flatfield correction, see also FlatFieldCorrection_20170320.m
%   Detailed explanation goes here

    %vol(CorrIm2D)
    for i = 1:size(RawIm3D, 3)
        RawIm3D(:,:,i) = double(RawIm3D(:,:,i)) .* double(CorrIm2D);
    end
    %vol(RawIm3D,0,1000,'hot')
    RawIm3D = uint16(RawIm3D);
    
end

