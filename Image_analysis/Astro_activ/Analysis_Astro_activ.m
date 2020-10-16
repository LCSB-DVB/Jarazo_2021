function [Objects] = Analysis_Astro_activ(InfoTableThis,PreviewPath, ch1, ch2, ch3, ch4)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

   %% image analysis

%% Segment nuclei
                  
                    NucChannelDoG = imfilter(ch1, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(ch1, fspecial('gaussian', 10, 6), 'symmetric'); %it(NucChannelDoG);
                     NucMask = NucChannelDoG > 15;
                     NucMask = bwareaopen(NucMask, 10); %it(NucMask); it(ch1);

%% Segment S100bns
    S100bSegmentationIm = ch2; % S100b  % vol(ch4,0,250) it(ch2);
    S100bSegmentationImGlobal = imfilter(S100bSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric'); % it(S100bSegmentationImGlobal)
    S100bSegmentationImLocal = imfilter(S100bSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(S100bSegmentationIm, fspecial('gaussian', 5, 6), 'symmetric'); % it(S100bSegmentationImLocal)
    S100bMaskGlobal = S100bSegmentationImGlobal > 50;  % it(S100bMaskGlobal)
    S100bMaskLocal = S100bSegmentationImLocal > 11;  % it(S100bMaskLocal)
    S100bMask = S100bMaskGlobal | S100bMaskLocal; % Fill gaps % vol(S100bMask) it(S100bMask)

    % Remove small objects which are certainly no S100bns
    S100bMask = bwareaopen(S100bMask, 5);  
    
%% Segment TH
    THSegmentationIm = ch3; % it(ch3)
    THSegmentationImGlobal = imfilter(THSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric'); % it(THSegmentationImGlobal)
    THSegmentationImLocal = imfilter(THSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(THSegmentationIm, fspecial('gaussian', 5, 6), 'symmetric'); %  it(THSegmentationImLocal)
    THMaskGlobal = THSegmentationImGlobal > 35; % it(THMaskGlobal)
    THMaskLocal = THSegmentationImLocal > 2; % it(THMaskLocal)
    THMask = THMaskGlobal | THMaskLocal; % Fill gaps % it(THMask)

    % Remove small objects which are certainly no TH neurons
    THMask = bwareaopen(THMask, 5); %it(THMask) it(ch3)
    
    %% Segment Astro
    GFAPSegmentationIm = ch4; % GFAP  % it(ch4);
    GFAPSegmentationImGlobal = imfilter(GFAPSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric'); % it(GFAPSegmentationImGlobal)
    GFAPSegmentationImLocal = imfilter(GFAPSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(GFAPSegmentationIm, fspecial('gaussian', 5, 6), 'symmetric'); % it(GFAPSegmentationImLocal)
    GFAPMaskGlobal = GFAPSegmentationImGlobal > 40;  % it(GFAPMaskGlobal)
    GFAPMaskLocal = GFAPSegmentationImLocal > 8;  % it(GFAPMaskLocal)
    GFAPMask = GFAPMaskGlobal | GFAPMaskLocal; % Fill gaps % it(GFAPMask)

    % Remove small objects which are certainly no Astro
    GFAPMask = bwareaopen(GFAPMask, 10);  
    
    %% Colocalization GFAP S100b
    ColoGFAPS100b=S100bMask & GFAPMask; %it(ColoGFAPS100b)
    
      %% Collect data
        Objects = table();        
        InfoTableThisSample = InfoTableThis; 
        THinTHMask = ch3 .* uint16(THMask);
        S100binS100bMask = ch2 .* uint16(S100bMask);
        GFAPinGFAPMask = ch4 .* uint16(GFAPMask);
        Objects.THbyNucVol = sum(THinTHMask(:)) / sum(NucMask(:));
        Objects.THbyTHVol = sum(THinTHMask(:)) / sum(THMask(:)); % Average pixel value
        Objects.S100bbyS100bVol = sum(S100binS100bMask(:)) / sum(S100bMask(:));% Average pixel value
        Objects.S100bbyNucVol = sum(S100binS100bMask(:)) / sum(NucMask(:));
        Objects.GFAPbyGFAPVol = sum(GFAPinGFAPMask(:)) / sum(GFAPMask(:));% Average pixel value
        Objects.GFAPbyNucVol = sum(GFAPinGFAPMask(:)) / sum(NucMask(:));
        Objects.THVolByS100bVol = sum(THMask(:)) / sum(S100bMask(:)); 
        Objects.S100bVolByNucVol = sum(S100bMask(:)) / sum(NucMask(:));
        Objects.THVolByNucVol = sum(THMask(:)) / sum(NucMask(:));
        Objects.GFAPVolByNucVol = sum(GFAPMask(:)) / sum(NucMask(:));
        Objects.GFAPVolByTHVol = sum(GFAPMask(:)) / sum(THMask(:));
        Objects.GFAPVolByS100bVol = sum(GFAPMask(:)) / sum(S100bMask(:));
        Objects.ColoGFAPS100bVolByNucVol = sum(ColoGFAPS100b(:)) / sum(NucMask(:));
        Objects.ColoGFAPS100bVolByS100bVol = sum(ColoGFAPS100b(:)) / sum(THMask(:)); 
        Objects.S100bVol = sum(S100bMask(:)); % pixel count
        Objects.THVol = sum(THMask(:));
        Objects.NucVol = sum(NucMask(:));
        Objects.GFAPVol = sum(GFAPMask(:));
        Objects.ColoGFAPS100bVol = sum(ColoGFAPS100b(:));
        Objects = [InfoTableThisSample, Objects];
        
    %% save 2D previews
    
    
    Barcode = InfoTableThis.Barcode{:};
    Areaname = InfoTableThis.AreaName{:}; 
    row = InfoTableThis.Row;
    column = InfoTableThis.Column;
    field = InfoTableThis.field;
    
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, ~] = f_barMask(50, 0.646, imSize, imSize(1)-50, 550, 10);
    %it(BarMask)
    RGB_1 = cat(3, imadjust(max(ch4, [], 3), [0 0.002], [0 1]), imadjust(max(ch2, [], 3), [0.0005 0.0015], [0 1]),imadjust(max(ch1, [], 3), [0 0.005], [0 1]));
    RGB_1 = imoverlay(RGB_1, BarMask, [1 1 1]); % imtool(RGB_1)
    
    Previewch4 = imadjust(max(ch4, [], 3), [0 0.002], [0 1]);
    Previewch1 = imadjust(max(ch1, [], 3), [0 0.005], [0 1]);
    Previewch3 = imadjust(max(ch3, [], 3), [0 0.002], [0 1]);
    Previewch2 = imadjust(max(ch2, [], 3), [0.0004 0.0015], [0 1]);
    %it(Previewch4)
    %it(Previewch1)
    %it(Previewch3)
    %it(Previewch2)
    
    % Channels with Mask
    PreviewTH = imoverlay2(Previewch3, bwperim(max(THMask,[],3)), [1 0 0]);
    %PreviewTH = imoverlay2(PreviewTH, BarMask, [1 1 1]);
    PreviewS100b = imoverlay2(Previewch2, bwperim(max(S100bMask,[],3)), [0 1 0]);
    %PreviewS100b = imoverlay2(PreviewS100b, BarMask, [1 1 1]);
    PreviewNuc = imoverlay2(Previewch1, bwperim(max(NucMask,[],3)), [0 0 1]);
    %PreviewNuc = imoverlay2(PreviewNuc, BarMask, [1 1 1]);
    PreviewGFAP = imoverlay2(Previewch4, bwperim(max(GFAPMask,[],3)), [1 0 0]);
    %PreviewGFAP = imoverlay2(PreviewGFAP, BarMask, [1 1 1]);
    %it(PreviewTH)
    %it(PreviewS100b)
    %it(PreviewNuc)
    %it(PreviewGFAP)
          
    filename_ch3_TH = sprintf('%s\\%s_%s_%03d%03d%03d_ch3_TH.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch1_Nuc = sprintf('%s\\%s_%s_%03d%03d%03d_ch1_Nuc.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch2_S100b= sprintf('%s\\%s_%s_%03d%03d%03d_ch2_S100b.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch4_GFAP= sprintf('%s\\%s_%s_%03d%03d%03d_ch4_GFAP.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_THMask = sprintf('%s\\%s_%s_%03d%03d%03d_THMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_S100bMask = sprintf('%s\\%s_%s_%03d%03d%03d_S100bMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_NucMask = sprintf('%s\\%s_%s_%03d%03d%03d_NucMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_GFAPMask = sprintf('%s\\%s_%s_%03d%03d%03d_GFAPMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_RGB_1 = sprintf('%s\\%s_%s_%03d%03d%03d_GFAP_S100b.png', PreviewPath, Barcode, Areaname, row, column, field);  
    
    imwrite(Previewch3, filename_ch3_TH);
    imwrite(Previewch1, filename_ch1_Nuc);
    imwrite(Previewch2, filename_ch2_S100b);
    imwrite(Previewch4, filename_ch4_GFAP);
    imwrite(PreviewS100b, filename_S100bMask);
    imwrite(PreviewNuc, filename_NucMask);
    imwrite(PreviewTH, filename_THMask);
    imwrite(PreviewGFAP, filename_GFAPMask);
    imwrite(RGB_1, filename_RGB_1);
   
    
end

