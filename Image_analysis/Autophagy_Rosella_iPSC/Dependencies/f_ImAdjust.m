function [ImAdjusted] = f_ImAdjust(Im, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %ImAdjusted = Im / max(Im(:)) * ((2^16)-1);
    if nargin == 1
        ImAdjusted = double(Im) / double(quantile(Im(:), 0.99)) * ((2^16)-1);
        ImAdjusted = cast(ImAdjusted, class(Im));
    elseif nargin == 2
        ImAdjusted = double(Im) / double(quantile(Im(:), varargin{1})) * ((2^16)-1);
        ImAdjusted = cast(ImAdjusted, class(Im));
    end
    
end

