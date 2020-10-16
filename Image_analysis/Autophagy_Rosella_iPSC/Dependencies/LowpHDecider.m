function [Decision] = LowpHDecider(x,y,z)
%Check if connected components awere recognized in the autophagosome mask
%and if they ar sufficiently circular
%   x: 'autophagosomeMaskValues' 
%   y: 'NonCircularMaskValues'

    AutophagosomePixels = max(x{:});
    NonCircularClass = max(y{:});
    RatioFirstMedian = median(z{:});

    %if   (AutophagosomePixels > 0) & (NonCircularClass == 0)
    if   (AutophagosomePixels > 0) & (NonCircularClass == 0) | (RatioFirstMedian < 2)
        Decision = 1;
    else
        Decision = 0;
    end

end
