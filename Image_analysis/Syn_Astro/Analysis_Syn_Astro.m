function [Objects] = Analysis_Syn_Astro_for_FP(InfoTableThis,PreviewPath, ch1, ch2, ch3, ch4)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

   %% image analysis

%% Segment nuclei
                  
                    NucChannelDoG = imfilter(ch1, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(ch1, fspecial('gaussian', 10, 6), 'symmetric'); %it(NucChannelDoG);
                     NucMask = NucChannelDoG > 15;
                     NucMask = bwareaopen(NucMask, 10); %it(NucMask); it(ch1);

%% Segment Synns
    SynSegmentationIm = ch2; % Syn  % vol(ch4,0,250) it(ch2);
    SynSegmentationImGlobal = imfilter(SynSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric'); % it(SynSegmentationImGlobal)
    SynSegmentationImLocal = imfilter(SynSegmentationIm, fspecial('gaussian', 1, 3), 'symmetric')...
        - imfilter(SynSegmentationIm, fspecial('gaussian', 5, 6), 'symmetric'); % it(SynSegmentationImLocal)
    SynMaskGlobal = SynSegmentationImGlobal > 120;  % it(SynMaskGlobal)
    SynMaskLocal = SynSegmentationImLocal > 20;  % it(SynMaskLocal)
    SynMask = SynMaskGlobal | SynMaskLocal; % Fill gaps % vol(SynMask) it(SynMask)

    % Remove small objects which are certainly no Synns
    SynMask = bwareaopen(SynMask, 3);  
    
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
    
    %% Colocalization GFAP Syn
    ColoTHSyn=SynMask & THMask; %it(ColoTHSyn)
    ColoGFAPSyn=SynMask & GFAPMask;
      
	  %% Collect data
        Objects = table();        
        InfoTableThisSample = InfoTableThis; 
        THinTHMask = ch3 .* uint16(THMask);
        SyninSynMask = ch2 .* uint16(SynMask);
        GFAPinGFAPMask = ch4 .* uint16(GFAPMask);
        Objects.THbyNucVol = sum(THinTHMask(:)) / sum(NucMask(:));
        Objects.THbyTHVol = sum(THinTHMask(:)) / sum(THMask(:)); % Average pixel value
        Objects.SynbySynVol = sum(SyninSynMask(:)) / sum(SynMask(:));% Average pixel value
        Objects.SynbyNucVol = sum(SyninSynMask(:)) / sum(NucMask(:));
        Objects.GFAPbyGFAPVol = sum(GFAPinGFAPMask(:)) / sum(GFAPMask(:));% Average pixel value
        Objects.GFAPbyNucVol = sum(GFAPinGFAPMask(:)) / sum(NucMask(:));
        Objects.THVolBySynVol = sum(THMask(:)) / sum(SynMask(:)); 
        Objects.SynVolByNucVol = sum(SynMask(:)) / sum(NucMask(:));
        Objects.THVolByNucVol = sum(THMask(:)) / sum(NucMask(:));
        Objects.GFAPVolByNucVol = sum(GFAPMask(:)) / sum(NucMask(:));
        Objects.GFAPVolByTHVol = sum(GFAPMask(:)) / sum(THMask(:));
        Objects.GFAPVolBySynVol = sum(GFAPMask(:)) / sum(SynMask(:));
        Objects.ColoTHSynVolByNucVol = sum(ColoTHSyn(:)) / sum(NucMask(:));
        Objects.ColoTHSynVolBySynVol = sum(ColoTHSyn(:)) / sum(THMask(:));
		Objects.ColoGFAPSynVolByNucVol = sum(ColoGFAPSyn(:)) / sum(NucMask(:));
        Objects.ColoGFAPSynVolBySynVol = sum(ColoGFAPSyn(:)) / sum(THMask(:));		
        Objects.SynVol = sum(SynMask(:)); % pixel count
        Objects.THVol = sum(THMask(:));
        Objects.NucVol = sum(NucMask(:));
        Objects.GFAPVol = sum(GFAPMask(:));
        Objects.ColoTHSynVol = sum(ColoTHSyn(:));
		Objects.ColoGFAPSynVol = sum(ColoGFAPSyn(:));
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
    RGB = cat(3, imadjust(max(ch3, [], 3), [0 0.002], [0 1]), imadjust(max(ch2, [], 3), [0.00045 0.0015], [0 1]),imadjust(max(ch1, [], 3), [0 0.005], [0 1]) );
    RGB = imoverlay(RGB, BarMask, [1 1 1]); % imtool(RGB)
       
    Previewch1 = imadjust(max(ch1, [], 3), [0 0.005], [0 1]);
    Previewch3 = imadjust(max(ch3, [], 3), [0 0.002], [0 1]);
    Previewch2 = imadjust(max(ch2, [], 3), [0.00045 0.002], [0 1]);
    %it(Previewch4)
    %it(Previewch1)
    %it(Previewch3)
    %it(Previewch2)
    
    
    % Channels with Mask
    PreviewTH = imoverlay2(Previewch3, bwperim(max(THMask,[],3)), [1 0 0]);
    PreviewSyn = imoverlay2(Previewch2, bwperim(max(SynMask,[],3)), [0 1 0]);
    PreviewNuc = imoverlay2(Previewch1, bwperim(max(NucMask,[],3)), [0 0 1]);
    %it(PreviewTH)
    %it(PreviewSyn)
    %it(PreviewNuc)
           
    filename_ch3_TH = sprintf('%s\\%s_%s_%03d%03d%03d_ch3_TH.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch1_Nuc = sprintf('%s\\%s_%s_%03d%03d%03d_ch1_Nuc.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_ch2_Syn= sprintf('%s\\%s_%s_%03d%03d%03d_ch2_Syn.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_THMask = sprintf('%s\\%s_%s_%03d%03d%03d_THMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_SynMask = sprintf('%s\\%s_%s_%03d%03d%03d_SynMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_NucMask = sprintf('%s\\%s_%s_%03d%03d%03d_NucMask.png', PreviewPath, Barcode, Areaname, row, column, field);
    filename_RGB = sprintf('%s\\%s_%s_%03d%03d%03d_TH_Syn.png', PreviewPath, Barcode, Areaname, row, column, field);
       
    imwrite(Previewch3, filename_ch3_TH);
    imwrite(Previewch1, filename_ch1_Nuc);
    imwrite(Previewch2, filename_ch2_Syn);
    imwrite(PreviewSyn, filename_SynMask);
    imwrite(PreviewNuc, filename_NucMask);
    imwrite(PreviewTH, filename_THMask);
    imwrite(RGB, filename_RGB);
 
   if sum(ColoTHSyn(:))> 0
      RGB_ColoTHSyn = cat(3, imadjust(max(ch3, [], 3), [0 0.002], [0 1]), imadjust(max(ch2, [], 3), [0.00045 0.002], [0 1]),imadjust(max(ch1, [], 3), [0 0.005], [0 1]) );
      RGB_ColoTHSyn = imoverlay2(RGB_ColoTHSyn, max(ColoTHSyn,[],3), [1 1 1]); 
      PreviewColoTHPath = [PreviewPath, filesep, 'ColoTH'];
      filename_RGB_ColoTHSyn = sprintf('%s\\%s_%s_%03d%03d%03d_RGB_ColoTH.png', PreviewColoTHPath, Barcode, Areaname, row, column, field);
      imwrite(RGB_ColoTHSyn, filename_RGB_ColoTHSyn);
    else
   end
       
end

