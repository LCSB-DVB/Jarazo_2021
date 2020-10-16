function [Objects] = ImageAnalysis_S3_auto_modul(InfoTableThis,DataPath, ch1, ch2, ch3, ch4)

   %% image analysis

%% Segment nuclei
                  
                    NucChannelDoG = imfilter(ch1, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(ch1, fspecial('gaussian', 10, 6), 'symmetric'); %it(NucChannelDoG);
                     NucMask = NucChannelDoG > 15;
                     NucMask = bwareaopen(NucMask, 10); %it(NucMask); it(ch1);

%% Segment neurons
    NeuroSegmentationIm = ch2; % Tuj1  % vol(ch4,0,250) it(ch2);
    NeuroSegmentationImGlobal = imfilter(NeuroSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric'); % it(NeuroSegmentationImGlobal)
    NeuroSegmentationImLocal = imfilter(NeuroSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(NeuroSegmentationIm, fspecial('gaussian', 5, 6), 'symmetric'); % it(NeuroSegmentationImLocal)
    NeuroMaskGlobal = NeuroSegmentationImGlobal > 700;  % it(NeuroMaskGlobal)
    NeuroMaskLocal = NeuroSegmentationImLocal > 100;  % it(NeuroMaskLocal)
    NeuroMask = NeuroMaskGlobal | NeuroMaskLocal; % Fill gaps % vol(NeuroMask) it(NeuroMask)

    % Remove small objects which are certainly no neurons
    NeuroMask = bwareaopen(NeuroMask, 5);  
    
%% Segment TH
    THSegmentationIm = ch3; % it(ch3)
    THSegmentationImGlobal = imfilter(THSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric'); % it(THSegmentationImGlobal)
    THSegmentationImLocal = imfilter(THSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(THSegmentationIm, fspecial('gaussian', 5, 6), 'symmetric'); %  it(THSegmentationImLocal)
    THMaskGlobal = THSegmentationImGlobal > 50; % it(THMaskGlobal)
    THMaskLocal = THSegmentationImLocal > 8; % it(THMaskLocal)
    THMask = THMaskGlobal | THMaskLocal; % Fill gaps % it(THMask)

    % Remove small objects which are certainly no TH neurons
    THMask = bwareaopen(THMask, 3); %it(THMask) it(ch3)
    
    %% Segment Astro
    GFAPSegmentationIm = ch4; % GFAP  % it(ch4);
    GFAPSegmentationImGlobal = imfilter(GFAPSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric'); % it(GFAPSegmentationImGlobal)
    GFAPSegmentationImLocal = imfilter(GFAPSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(GFAPSegmentationIm, fspecial('gaussian', 5, 6), 'symmetric'); % it(GFAPSegmentationImLocal)
    GFAPMaskGlobal = GFAPSegmentationImGlobal > 400;  % it(GFAPMaskGlobal)
    GFAPMaskLocal = GFAPSegmentationImLocal > 30;  % it(GFAPMaskLocal)
    GFAPMask = ~NucMask & (GFAPMaskGlobal | GFAPMaskLocal) ; % Fill gaps % it(GFAPMask)

    % Remove small objects which are certainly no Astro
    GFAPMask = bwareaopen(GFAPMask, 10);  
    
    %% Colocalization TH Tuj1
    ColoTHTuj1=NeuroMask & THMask; %it(ColoTHTuj1)
    
      %% Collect data
        Objects = table();        
        InfoTableThisSample = InfoTableThis; 
        THinTHMask = ch3 .* uint16(THMask);
        Tuj1inTuj1Mask = ch2 .* uint16(NeuroMask);
        GFAPinGFAPMask = ch4 .* uint16(GFAPMask);
        Objects.THbyNucVol = sum(THinTHMask(:)) / sum(NucMask(:));
        Objects.THbyTHVol = sum(THinTHMask(:)) / sum(THMask(:)); % Average pixel value
        Objects.Tuj1byTuj1Vol = sum(Tuj1inTuj1Mask(:)) / sum(NeuroMask(:));% Average pixel value
        Objects.Tuj1byNucVol = sum(Tuj1inTuj1Mask(:)) / sum(NucMask(:));
        Objects.GFAPbyGFAPVol = sum(GFAPinGFAPMask(:)) / sum(GFAPMask(:));% Average pixel value
        Objects.GFAPbyNucVol = sum(GFAPinGFAPMask(:)) / sum(NucMask(:));
        Objects.THVolByTuj1Vol = sum(THMask(:)) / sum(NeuroMask(:)); 
        Objects.Tuj1VolByNucVol = sum(NeuroMask(:)) / sum(NucMask(:));
        Objects.THVolByNucVol = sum(THMask(:)) / sum(NucMask(:));
        Objects.GFAPVolByNucVol = sum(GFAPMask(:)) / sum(NucMask(:));
        Objects.GFAPVolByTHVol = sum(GFAPMask(:)) / sum(THMask(:));
        Objects.GFAPVolByTuj1Vol = sum(GFAPMask(:)) / sum(NeuroMask(:));
        Objects.ColoTHTuj1VolByNucVol = sum(ColoTHTuj1(:)) / sum(NucMask(:));
        Objects.ColoTHTuj1VolByTuj1Vol = sum(ColoTHTuj1(:)) / sum(NeuroMask(:)); 
        Objects.Tuj1Vol = sum(NeuroMask(:)); % pixel count
        Objects.THVol = sum(THMask(:));
        Objects.NucVol = sum(NucMask(:));
        Objects.GFAPVol = sum(GFAPMask(:));
        Objects.ColoTHTuj1Vol = sum(ColoTHTuj1(:));
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
    [BarMask, ~] = f_barMask(50, 0.646, imSize, imSize(1)-50, 550, 10);
    %it(BarMask)
    RGB = cat(3, imadjust(max(ch3, [], 3)), imadjust(max(ch2, [], 3)),imadjust(max(ch1, [], 3)) );
    RGB = imoverlay(RGB, BarMask, [1 1 1]); % imtool(RGB)
      
    Previewch4 = imadjust(max(ch4,[],3));
    Previewch1 = imadjust(max(ch1,[],3));
    Previewch3 = imadjust(max(ch3,[],3));
    Previewch2 = imadjust(max(ch2,[],3));
    %it(Previewch4)
    %it(Previewch1)
    %it(Previewch3)
    
    % Channels with Mask
    PreviewTH = imoverlay2(imadjust(max(ch3,[],3)), bwperim(max(THMask,[],3)), [1 0 0]);
    PreviewTH = imoverlay2(PreviewTH, BarMask, [1 1 1]);
    PreviewTuj1 = imoverlay2(imadjust(max(ch2,[],3)), bwperim(max(NeuroMask,[],3)), [0 1 0]);
    PreviewTuj1 = imoverlay2(PreviewTuj1, BarMask, [1 1 1]);
    PreviewNuc = imoverlay2(imadjust(max(ch1,[],3)), bwperim(max(NucMask,[],3)), [0 0 1]);
    PreviewNuc = imoverlay2(PreviewNuc, BarMask, [1 1 1]);
    PreviewGFAP = imoverlay2(imadjust(max(ch4,[],3)), bwperim(max(GFAPMask,[],3)), [1 0 0]);
    PreviewGFAP = imoverlay2(PreviewGFAP, BarMask, [1 1 1]);
    %it(PreviewTH)
    %it(PreviewTuj1)
    %it(PreviewNuc)
    %it(PreviewGFAP)
          
    filename_ch3_TH = sprintf('%s\\%s_%s_%03d%03d%03d_ch3_TH.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch1_Nuc = sprintf('%s\\%s_%s_%03d%03d%03d_ch1_Nuc.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch2_Tuj1= sprintf('%s\\%s_%s_%03d%03d%03d_ch2_Tuj1.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch4_GFAP= sprintf('%s\\%s_%s_%03d%03d%03d_ch4_GFAP.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_THMask = sprintf('%s\\%s_%s_%03d%03d%03d_THMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_Tuj1Mask = sprintf('%s\\%s_%s_%03d%03d%03d_Tuj1Mask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_NucMask = sprintf('%s\\%s_%s_%03d%03d%03d_NucMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_GFAPMask = sprintf('%s\\%s_%s_%03d%03d%03d_GFAPMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_RGB = sprintf('%s\\%s_%s_%03d%03d%03d_TH_Tuj1.png', PreviewPath, Barcode, Areaname, row, column, field);
    
    
    imwrite(Previewch3, filename_ch3_TH);
    imwrite(Previewch1, filename_ch1_Nuc);
    imwrite(Previewch2, filename_ch2_Tuj1);
    imwrite(Previewch4, filename_ch4_GFAP);
    imwrite(PreviewTuj1, filename_Tuj1Mask);
    imwrite(PreviewNuc, filename_NucMask);
    imwrite(PreviewTH, filename_THMask);
    imwrite(PreviewGFAP, filename_GFAPMask);
    imwrite(RGB, filename_RGB);
   
    
end

