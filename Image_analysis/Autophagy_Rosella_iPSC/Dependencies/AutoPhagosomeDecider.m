function [Decision] = AutoPhagosomeDecider(x,y)
%Check if connected components awere recognized in the autophagosome mask
%and if they ar sufficiently circular
%   x: 'autophagosomeMaskValues' 
%   y: 'NonCircularMaskValues'

    AutophagosomePixels = max(x{:});
    NonCircularClass = max(y{:});

    if   (AutophagosomePixels > 0) & (NonCircularClass == 0)
        Decision = 1;
    else
        Decision = 0;
    end

end

