function [Decision] = LateAutoysosomeDecider(x,y,z,l,omega,CenterGreen,OutsideGreen)
%After selecting Autophagosomes and Phagophores, the function is checking 
%both Masks in Green and Red and deciding based on proportion to classify
%Late Autolysosomes
%   x: 'redMaskValues' 
%   y: 'greenMaskValues'
%   z: 'AutophagosomeDecision'
%   l: 'PhagophoreDecision'



    PixelCount = numel(x{:});
    % Check red signal
    RedProportion = sum(x{:})/PixelCount;

    % Check green signal
    GreenProportion = sum(y{:})/PixelCount;
    
    % Check ratio image
    RatioVec = omega{:};
    MedianRatio = median(RatioVec);

    if ((RedProportion > 0.25) & (GreenProportion < 0.1) & (z == 0) & (l == 0)) | ((MedianRatio <= 2) & (z == 0) & (l == 0))
        Decision = 1;
%     elseif CenterGreen >= PerimGreen % Early autolysosomes in contrast to late autolysosomes have a bright green ring
%         Decision = 1;
    elseif CenterGreen <= OutsideGreen % Early autolysosomes in contrast to late autolysosomes have a bright green ring
        Decision = 1;

    else
        Decision = 0;
    end

end

