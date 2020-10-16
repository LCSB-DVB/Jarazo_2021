function [Objects] = Analysis_Diff(InfoTableThis,DataPath, ch1, ch2, ch3, ch4)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% image analysis

%% Segment nuclei
                  
                    NucChannelDoG = imfilter(ch1, fspecial('gaussian', 10, 2), 'symmetric')...
        - imfilter(ch1, fspecial('gaussian', 20, 6), 'symmetric'); %it(NucChannelDoG);
                     NucMask = NucChannelDoG > 15;
                     NucMask = bwareaopen(NucMask, 10); %it(NucMask); it(ch1);

%% Segment neurons
    NeuroSegmentationIm = ch4; % Tuj1 only %it(ch4);
          
    NeuroSegmentationImGlobal = imfilter(NeuroSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric'); % imtool(NeuroSegmentationImGlobal)
    NeuroSegmentationImLocal = imfilter(NeuroSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(NeuroSegmentationIm, fspecial('gaussian', 5, 6), 'symmetric'); % vol(NeuroSegmentationImLocal,0,20)
 
    NeuroMaskGlobal = NeuroSegmentationImGlobal > 100; % vol(NeuroSegmentationImGlobal) % vol(NeuroMaskGlobal)
   
    NeuroMaskLocal = NeuroSegmentationImLocal > 2; % vol(NeuroSegmentationImLocal) % vol(NeuroMaskLocal)
    
    NeuroMask = NeuroMaskGlobal | NeuroMaskLocal; % Fill gaps % vol(NeuroMask) it(NeuroMask)

    % Remove small objects which are certainly no neurons
    NeuroMask = bwareaopen(NeuroMask, 5);                        
%% Segment TH
    THSegmentationIm = ch3; % it(ch3)
    %vol(THSegmentationIm, 0,50)
    THSegmentationImGlobal = imfilter(THSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric'); % it(THSegmentationImGlobal)
    THSegmentationImLocal = imfilter(THSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(THSegmentationIm, fspecial('gaussian', 5, 6), 'symmetric'); % vol(THSegmentationImLocal,0,20) it(THSegmentationImLocal)
    THMaskGlobal = THSegmentationImGlobal > 100; % vol(THSegmentationImGlobal) % vol(THMaskGlobal) it(THMaskGlobal)
    THMaskLocal = THSegmentationImLocal > 5; % vol(THSegmentationImLocal) % vol(THMaskLocal) it(THMaskLocal)
    THMask = THMaskGlobal | THMaskLocal; % Fill gaps % vol(THMask)

    % Remove small objects which are certainly no TH neurons
    THMask = bwareaopen(THMask, 10); %it(THMask) it(ch3)

    
      %% Collect data
        Objects = table();        
        InfoTableThisSample = InfoTableThis; 
        THinTHMask = ch1 .* uint16(THMask);
        Tuj1inTuj1Mask = ch3 .* uint16(NeuroMask);
        Objects.THbyNucVol = sum(THinTHMask(:)) / sum(NucMask(:));
        Objects.THbyTHVol = sum(THinTHMask(:)) / sum(THMask(:)); % Average pixel value
        Objects.Tuj1byTuj1Vol = sum(Tuj1inTuj1Mask(:)) / sum(NeuroMask(:));
        Objects.Tuj1byNucVol = sum(Tuj1inTuj1Mask(:)) / sum(NucMask(:));
        Objects.THVolByTuj1Vol = sum(THMask(:)) / sum(NeuroMask(:)); % Average pixel value
        Objects.Tuj1VolByNucVol = sum(NeuroMask(:)) / sum(NucMask(:));
        Objects.THVolByNucVol = sum(THMask(:)) / sum(NucMask(:));
        Objects.Tuj1Vol = sum(NeuroMask(:)); % pixel count
        Objects.THVol = sum(THMask(:));
        Objects.NucVol = sum(NucMask(:));
        Objects = [InfoTableThisSample, Objects];
        
         %% save 2D previews
    
    PreviewPath = [DataPath, filesep, 'Previews'];
    Barcode = InfoTableThis.Barcode{:};
    Areaname = InfoTableThis.AreaName{:}; 
    row = InfoTableThis.Row;
    column = InfoTableThis.Column;
    field = InfoTableThis.field;
    
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask,~] = f_barMask(50, 0.646, imSize, imSize(1)-50, 550, 10);
    %it(BarMask)
    %RGB
    RGB = cat(3, imadjust(max(ch3, [], 3), [0 0.007], [0 1]), imadjust(max(ch4, [], 3), [0 0.0025], [0 1]),imadjust(max(ch1, [], 3), [0 0.009], [0 1]) );
    RGB = imoverlay(RGB, BarMask, [1 1 1]); % imtool(RGB)

    Previewch4 = imadjust(max(ch4,[],3));
    Previewch1 = imadjust(max(ch1,[],3));
    Previewch3 = imadjust(max(ch3,[],3));
    %it(Previewch4)
    %it(Previewch1)
    %it(Previewch3)
    
    % Channels with Mask
    PreviewTH = imoverlay2(imadjust(max(ch3,[],3)), bwperim(max(THMask,[],3)), [1 0 0]);
    PreviewTH = imoverlay2(PreviewTH, BarMask, [1 1 1]);
    PreviewTuj1 = imoverlay2(imadjust(max(ch4,[],3)), bwperim(max(NeuroMask,[],3)), [0 1 0]);
    PreviewTuj1 = imoverlay2(PreviewTuj1, BarMask, [1 1 1]);
    PreviewNuc = imoverlay2(imadjust(max(ch1,[],3)), bwperim(max(NucMask,[],3)), [0 0 1]);
    PreviewNuc = imoverlay2(PreviewNuc, BarMask, [1 1 1]);
    %it(PreviewTH)
    %it(PreviewTuj1)
    %it(PreviewNuc)
   
       
    filename_ch3_TH = sprintf('%s\\%s_%s_%03d%03d%03d_ch3_TH.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch1_Nuc = sprintf('%s\\%s_%s_%03d%03d%03d_ch1_Nuc.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch4_Tuj1= sprintf('%s\\%s_%s_%03d%03d%03d_ch4_Tuj1.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_THMask = sprintf('%s\\%s_%s_%03d%03d%03d_THMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_Tuj1Mask = sprintf('%s\\%s_%s_%03d%03d%03d_Tuj1Mask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_NucMask = sprintf('%s\\%s_%s_%03d%03d%03d_NucMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_RGB = sprintf('%s\\%s_%s_%03d%03d%03d_RGB.png', PreviewPath, Barcode, Areaname, row, column, field);
      
    
    imwrite(Previewch3, filename_ch3_TH);
    imwrite(Previewch1, filename_ch1_Nuc);
    imwrite(Previewch4, filename_ch4_Tuj1);
    imwrite(PreviewTuj1, filename_Tuj1Mask);
    imwrite(PreviewNuc, filename_NucMask);
    imwrite(PreviewTH, filename_THMask);
    imwrite(RGB, filename_RGB);
   
    
end

