function [ FigureHandle ] = f_imoverlay3D_5Masks(Im3D, Mask3D_R, Mask3D_G, Mask3D_B, Mask3D_V, Mask3D_C, Title)
%Create figure showing masks in color on 3D gray image
%   Detailed explanation goes here

% Normalize Im3D
Im3D = double(Im3D);
Im3D = Im3D ./ max(Im3D(:));

subplotsNumel = size(Im3D,3);
SideLength = ceil(sqrt(subplotsNumel));

%axCell = {};

FigureHandle = figure
for s = 1:subplotsNumel
    PlaneOverlay = imoverlay2(imadjust(Im3D(:,:,s)), Mask3D_R(:,:,s), [1 0 0]);
    PlaneOverlay = imoverlay2( PlaneOverlay, Mask3D_G(:,:,s), [0 1 0]);
    PlaneOverlay = imoverlay2( PlaneOverlay, Mask3D_B(:,:,s), [0 0 1]);
    PlaneOverlay = imoverlay2( PlaneOverlay, Mask3D_V(:,:,s), [1 0 1]);
    PlaneOverlay = imoverlay2( PlaneOverlay, Mask3D_C(:,:,s), [0 1 1]);
    ax(s) = subplot(SideLength,SideLength,s)
    imshow(PlaneOverlay), title(['Plane ' num2str(s)])
    %imshow(Im3D(:,:,s),[])
end

annotation('textbox', [0, 0.9, 1, 0.1], 'String', Title, 'EdgeColor', 'none', 'HorizontalAlignment', 'Center', 'FontSize', 14)
linkaxes(ax, 'xy');