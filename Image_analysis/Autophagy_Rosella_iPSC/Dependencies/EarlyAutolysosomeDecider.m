function [Decision] = EarlyAutolysosomeDecider(x,y,z,l,m,AutoLysoMask,CenterGreen,PerimGreen)
%After selecting Autophagosomes, Phagophores and Late Autolysosomes, the function is checking 
%both Masks in Green and Red and deciding based on proportion to classify
%Early Autolysosomes
%   x: 'redMaskValues' 
%   y: 'greenMaskValues'
%   z: 'AutophagosomeDecision'
%   l: 'PhagophoreDecision'
%   m: 'LateAutolysosomeDecision'


    PixelCount = numel(x{:});
    % Check red signal
    RedProportion = sum(x{:})/PixelCount;
    % AutoLysoCheck
    AutoLysoCheck = max(AutoLysoMask{:});

    % Check green signal
    GreenProportion = sum(y{:})/PixelCount;
    if ((RedProportion > 0.25) & (GreenProportion <= 0.25) & (z == 0) & (l == 0) & (m == 0)) | AutoLysoCheck
        Decision = 1;
    elseif PerimGreen < CenterGreen % the perimter has to be brighter than the center
        Decision = 0; % send to trash decision
    else
        Decision = 0;
    end

end

