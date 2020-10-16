function [MaskCleared] = f_imclearborder3D(Mask)
%Clear border along xy but not along z

    MaskCleared = padarray(Mask, [0 0 1], 0, 'both');
    MaskCleared = imclearborder(MaskCleared);
    MaskCleared = MaskCleared(:,:,2:end-1);

end

