function [SlidePreviewIms] = f_CV8000_SlideViewOrganoids(Path, WellsToDisplay, channel, GrayRangeInput)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    SlidePreviewIms = {};
    filesAll = sort(dirrec(Path, '.mat')');
    %Wells = regexp(filesAll, '.*\\(...)_organoidID.*', 'tokens');
    Wells = regexp(filesAll, '.*/(...)_organoidID.*', 'tokens');    
    Wells = cellfun(@(x) x{:}{:}, Wells, 'UniformOutput', false);
    %WellsSelect = strcmp(WellsSelect, Wells);
    %Wells = ismember(Wells, WellsToDisplay);
    %Wells = Wells;
    Slides = unique(Wells);
    ChannelsSelect = regexp(filesAll, '.*channel(.)', 'tokens');
    ChannelsSelect = cellfun(@(x) x{:}{:}, ChannelsSelect, 'UniformOutput', false);
    ChannelsSelect = strcmp(ChannelsSelect, num2str(channel));
    
    for Slide = 1:size(Slides, 1) % Here a slide corresponds to a well
        WellThis = Slides{Slide};
        SelectFiles = filesAll(strcmp(Wells, WellThis) & ChannelsSelect);
    
        %% Initialize preview
        FirstOrganoidData = load(SelectFiles{1});
        LabelMatrix = FirstOrganoidData.LabelMatrix; % imtool(LabelMatrix, [])
        OrganoidGroupIm = FirstOrganoidData.OrganoidGroupIm;
        SlidePreviewIm = uint8(padarray(LabelMatrix .* OrganoidGroupIm, [11, 11], 0, 'post'));
        SlidePreviewIm = imresize(SlidePreviewIm, 100/9, 'nearest');
        SlidePreviewIm = cat(3, SlidePreviewIm, SlidePreviewIm, SlidePreviewIm);%RGB
        SlidePreviewPatternIm = SlidePreviewIm;
        %imtool(SlidePreviewPatternIm,[])
        for i = 1:max(LabelMatrix(:))
            [r, c] = find(SlidePreviewPatternIm == i);
            CornerRow = min(r);
            CornerCol = min(c);
            OrganoidImThis = load(SelectFiles{i});
            OrganoidImThis = uint16(OrganoidImThis.OrganoidImThis);%imtool(OrganoidImThis)
            OrganoidImThisContour = imdilate(bwperim(medfilt2(OrganoidImThis) > 0), strel('disk', 2)); % imtool(imadjust(OrganoidImThis, [0 0.03], [0 1]))
            OrganoidImThis = imoverlay(imadjust(OrganoidImThis, GrayRangeInput, [0 1]), OrganoidImThisContour, [0 1 0]);%imtool(OrganoidImThis)
            OrganoidImThis = insertText(OrganoidImThis, [1 1], num2str(i), 'FontSize', 30);
            SlidePreviewIm(CornerRow:CornerRow+size(OrganoidImThis, 1)-1, CornerCol:CornerCol+size(OrganoidImThis, 2)-1,:) = OrganoidImThis;
        end

        SlidePreviewIm = insertText(SlidePreviewIm, [1 1], WellThis, 'FontSize', 60, 'BoxColor', 'red', 'BoxOpacity', 0.3, 'TextColor','white');
        %imtool(SlidePreviewIm,[])
        SlidePreviewIms{Slide, 1} = SlidePreviewIm;

    end
end

