function [Decision] = PhagophoreDecider(x,y,z,omega,ImG,ImR,AutoLysoMask, CenterGreen, OutsideGreen)
%After selecting Autophagosomes, the function is checking 
%both Masks in Green and Red and deciding based on proportion to classify
%Phagophores
%   x: 'redMaskValues' 
%   y: 'greenMaskValues'
%   z: 'AutophagosomeDecision'

%test

    PixelCount = numel(x{:});
    % Check red signal
    RedProportion = sum(x{:})/PixelCount;

    % Check green signal
    GreenProportion = sum(y{:})/PixelCount;
    
    % Check ratio image
    RatioVec = omega{:};
    MedianRatio = median(RatioVec);
    
    % Check Green
    Green = quantile(double(ImG{:}),0.75);
    
    % Check Red
    Red = quantile(double(ImR{:}),0.75);
    
    % AutoLyso check
    AutoLysoLabel = max(AutoLysoMask{:});
    
    

    %if   (GreenProportion > 0.25)&(RedProportion > 0.25) & (z == 0)
    %if   ((GreenProportion > 0.25)&(RedProportion > 0.25) & (z == 0)) & (~(MedianRatio < 2))
    if   ((GreenProportion > 0.25)&(RedProportion > 0.25) & (z == 0)) & (~(MedianRatio < 2)) & ~AutoLysoLabel
        Decision = 1;
    %elseif (Green > 7500) & (Red > 7500)
    %elseif (Green > 7500) & (Red > 4000)  & (z == 0) %this would mean red is intense in phagophores
    elseif (Green > 7500) & (Red > 4000)  & (z == 0)  & ~AutoLysoLabel %this would mean red is intense in phagophores
    %elseif (Green > 7500) & (z == 0)
    %elseif (Green > 0) & (Red > 0)
        Decision = 1;
    elseif double(CenterGreen) > (1.25 * (double(OutsideGreen)))
        Decision = 1;
    else
        Decision = 0;
    end


end

