function [Objects] = Analysis_S3_CD_untag(ch1,ch2,ch3,ch4,InfoTableThis, PreviewPath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%it(ch1)
%it(ch2, 0, 500)
%it(ch3, 0, 500)
  
%% Segment nuclei
    NucleiBlurred = imfilter(ch1, fspecial('gaussian',5, 4)); % it(NucleiBlurred)
    NucleiMask = NucleiBlurred > 200; % it(NucleiMask)
        if sum(NucleiMask(:))== 0
            Objects = {}; 
            return
        end
    NucleiDeadMask = NucleiBlurred > 1000; % it(NucleiDeadMask)
    %% Segment p62somes
    p62DoG = imfilter(ch3, fspecial('gaussian', 1, 3), 'symmetric') - imfilter(ch2, fspecial('gaussian', 5, 6), 'symmetric'); % it(p62DoG)
    p62MaskLocal = p62DoG > 50; % it(p62MaskLocal)
    p62Mask = p62MaskLocal; % & p62GlobalMask;
    p62Mask = bwareaopen (p62Mask, 2); % it(p62Mask)
    %% Segment LAMPsomes
    LAMPDoG = imfilter(ch2, fspecial('gaussian', 1, 3), 'symmetric') - imfilter(ch3, fspecial('gaussian', 5, 6), 'symmetric'); % it(LAMPDoG)
    LAMPMaskLocal = LAMPDoG > 100; % it(LAMPMaskLocal)
    LAMPMask = LAMPMaskLocal; % & LAMPGlobalMask;
    LAMPMask = bwareaopen (LAMPMask, 2); % it(LAMPMask)

    %% Segment TH Soma
    THBlurred = imfilter(ch4, fspecial('gaussian',20, 4)); % it(THBlurred)
    THSomaMask = THBlurred > 50;  % it(THSomaMask)
    THSomaMask = bwareaopen(THSomaMask, 600);
	
	%% Segment TH Neurite
	THNeuriteDoG = imfilter(ch4, fspecial('gaussian', 7, 1), 'symmetric') - imfilter(ch4, fspecial('gaussian', 7, 2), 'symmetric');  %it(THNeuriteDoG)
    THNeuriteMask = THNeuriteDoG > 5;  % it(THNeuriteMask)
    THNeuriteMask = bwareaopen(THNeuriteMask, 7);
    THNeuriteMask = THNeuriteMask & ~THSomaMask;
	
	%% p62-LAMP Colocalization Mask
    ColoLpMask = p62Mask & LAMPMask; %it(ColoLpMask) 
	    
	%% p62-LAMP-THSoma Colocalization Mask
    ColoLpTSMask = p62Mask & LAMPMask & THSomaMask; %it(ColoLpTSMask) 
    
	%% p62-LAMP-THNeurite Colocalization Mask
    ColoLpTNMask = p62Mask & LAMPMask & THNeuriteMask; %it(ColoLpTNMask) 
	
	%% p62-THSoma Colocalization Mask
    ColopTSMask = p62Mask & THSomaMask; %it(ColopTSMask) 

	%% p62-THSoma Colocalization Mask
    ColopTNMask = p62Mask & THNeuriteMask; %it(ColopTNMask) 

	%% LAMP-THSoma Colocalization Mask
    ColoLATSMask = LAMPMask & THSomaMask; %it(ColoLATSMask) 
	
	%% LAMP-THNeurite Colocalization Mask
    ColoLATNMask = LAMPMask & THNeuriteMask; %it(ColoLATNMask) 
	
	%% Morphometrics
    p62LM = bwlabeln(p62Mask);
    p62Objects = regionprops('table', p62LM, ch3,{'Area','MajorAxisLength','MinorAxisLength','MeanIntensity'});
    LAMPLM = bwlabeln(LAMPMask);
    LAMPObjects = regionprops('table', LAMPLM, ch2,{'Area','MajorAxisLength','MinorAxisLength','MeanIntensity'});
    ColoLpLM = bwlabeln(ColoLpMask);
    ColoLpObjects = regionprops('table', ColoLpLM,{'Area','MajorAxisLength','MinorAxisLength'});
    ColoLpTSLM = bwlabeln(ColoLpTSMask);
    ColoLpTSObjects = regionprops('table', ColoLpTSLM, {'Area','MajorAxisLength','MinorAxisLength'});
	ColoLpTNLM = bwlabeln(ColoLpTNMask);
    ColoLpTNObjects = regionprops('table', ColoLpTNLM, {'Area','MajorAxisLength','MinorAxisLength'});
    ColopTSLM = bwlabeln(ColopTSMask);
    ColopTSObjects = regionprops('table', ColopTSLM, {'Area','MajorAxisLength','MinorAxisLength'});
	ColopTNLM = bwlabeln(ColopTNMask);
    ColopTNObjects = regionprops('table', ColopTNLM, {'Area','MajorAxisLength','MinorAxisLength'});
	ColoLATSLM = bwlabeln(ColoLATSMask);
    ColoLATSObjects = regionprops('table', ColoLATSLM, {'Area','MajorAxisLength','MinorAxisLength'});
	ColoLATNLM = bwlabeln(ColoLATNMask);
    ColoLATNObjects = regionprops('table', ColoLATNLM, {'Area','MajorAxisLength','MinorAxisLength'});
	
	%% Extract features
    Objects = table();
    Objects.File = InfoTableThis.files;
    Objects.Barcode = InfoTableThis.Barcode;
    Objects.AreaName = InfoTableThis.AreaName;
    Objects.ROW = InfoTableThis.Row;
    Objects.COL = InfoTableThis.Column;
    Objects.Field = InfoTableThis.field;
     %% Nucleus derived
    Objects.NucArea = sum(NucleiMask(:));
    Objects.NucDeadArea = sum(NucleiDeadMask(:));
	%% TH derived
    Objects.THSomaArea = sum(THSomaMask(:));
	Objects.THNeuriteArea = sum(THNeuriteMask(:));
    
    %% 2D previews
    voxelSizeX = 0.32;%40x Bin2
    % Scalebar
    imSize = size(ch1);
    [BarMask, ~] = f_barMask(20, voxelSizeX, imSize, imSize(1)-50, 75, 10);
    %it(BarMask)
    RGB = cat(3, imadjust(ch2), imadjust(ch3), imadjust(ch1));
    RGB = imoverlay(RGB, BarMask, [1 1 1]);
	RGB_1 = cat(3, imadjust(ch2), imadjust(ch3), imadjust(ch4));
    RGB_1 = imoverlay(RGB_1, BarMask, [1 1 1]);
    % imtool(RGB)
    p62MaskPreview = imoverlay(imadjust(max(ch3, [], 3)), bwperim(max(p62Mask, [], 3)), [1 0 0]);
    p62MaskPreview = imoverlay(p62MaskPreview, BarMask, [1 1 1]);
    %imtool(p62MaskPreview )
    p62RawPreview = imadjust(max(ch3, [], 3));
    p62RawPreview = imoverlay(p62RawPreview, BarMask, [1 1 1]);
    %imtool(p62RawPreview)
    
    LAMPMaskPreview = imoverlay(imadjust(max(ch2, [], 3)), bwperim(max(LAMPMask, [], 3)), [1 0 0]);
    LAMPMaskPreview = imoverlay(LAMPMaskPreview, BarMask, [1 1 1]);
    %imtool(LAMPMaskPreview )
    LAMPRawPreview = imadjust(max(ch2, [], 3));
    LAMPRawPreview = imoverlay(LAMPRawPreview, BarMask, [1 1 1]);
    %imtool(LAMPRawPreview)
    
   
    SavePathp62MaskPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_p62_mask.png'];
    SavePathp62rawPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_p62_raw.png'];
    SavePathLAMPMaskPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_LAMP_mask.png'];
    SavePathLAMPrawPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_LAMP_raw.png'];
    SavePathRGBPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RGB.png'];
    SavePathRGB_1Preview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RGB_1.png'];
    
   
    imwrite(p62MaskPreview, SavePathp62MaskPreview)
    imwrite(p62RawPreview,  SavePathp62rawPreview)
    imwrite(LAMPMaskPreview, SavePathLAMPMaskPreview)
    imwrite(LAMPRawPreview,  SavePathLAMPrawPreview)
    imwrite(RGB, SavePathRGBPreview)
    imwrite(RGB_1, SavePathRGB_1Preview)
    
    %% p62 derived
    if sum(p62Mask(:))> 0
        Objects.p62Area = sum(p62Mask(:));
        Objects.p62AreaNorm = sum(p62Mask(:))/sum(NucleiMask(:));
        voxelSizeX = 0.32;%Bin2
        voxelSizeY = 0.32; 
        Objects.Minp62Vol = min(p62Objects.Area) * voxelSizeX * voxelSizeY; 
        Objects.Maxp62Vol = max(p62Objects.Area) * voxelSizeX * voxelSizeY;
        Objects.Meanp62Vol = mean(p62Objects.Area) * voxelSizeX * voxelSizeY ;
        Objects.Stdp62Vol = std(p62Objects.Area) * voxelSizeX * voxelSizeY;
        Objects.Medp62Vol = median(p62Objects.Area) * voxelSizeX * voxelSizeY;
        Objects.Madp62Vol = mad(p62Objects.Area, 1) * voxelSizeX * voxelSizeY;
        Objects.Totp62Volume = sum(p62Mask(:)) * voxelSizeX * voxelSizeY;
        Objects.Totp62VolumeNorm = (sum(p62Mask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
        Objects.Countp62 = size(p62Objects, 1);
        Objects.Countp62Norm = size(p62Objects, 1)/sum(NucleiMask(:));                
        Objects.MedMajorAxisLenghtp62 = median(p62Objects.MajorAxisLength)* voxelSizeX;
        Objects.MedMinorAxisLenghtp62 = median(p62Objects.MinorAxisLength)* voxelSizeX;
        Objects.Medp62Intensity = median(p62Objects.MeanIntensity);
        % Shape
        Conn6Strel = {};
        Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
        Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
        Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
        Conn6Strel = logical(cat(3, Conn6Strel{:}));
        p62ErodedMask = imerode(p62Mask, Conn6Strel);
        p62PerimMask = (p62Mask - p62ErodedMask) > 0;

        % Additional feature images and vectors
        p62BodyLabelIm = bwlabeln(p62ErodedMask, 6);

        % Erosion derived
        Objects.p62PerimPixels = sum(p62PerimMask(:));
        Objects.p62BodyPixels = sum(p62ErodedMask(:)); 
        Objects.p62BodyCount = max(p62BodyLabelIm(:)); % Needed for invagination feature
        Objects.p62ShapeBySurface = Objects.p62BodyPixels / Objects.p62PerimPixels; % Roundness feature
        Objects.p62BodycountByp62count = Objects.p62BodyCount / Objects.Countp62; % Invagination feature
        else
        Objects.p62Area = 0;
        Objects.p62AreaNorm = 0;
        Objects.Minp62Vol = 0;
        Objects.Maxp62Vol = 0;
        Objects.Meanp62Vol = 0;
        Objects.Stdp62Vol = 0;
        Objects.Medp62Vol = 0;
        Objects.Madp62Vol = 0;
        Objects.Totp62Volume = 0;
        Objects.Totp62VolumeNorm = 0;
        Objects.Countp62 = 0;
        Objects.Countp62Norm = 0; 
        Objects.MedMajorAxisLenghtp62 = 0;
        Objects.MedMinorAxisLenghtp62 = 0;
        Objects.Medp62Intensity =0;
        % Erosion derived
        Objects.p62PerimPixels = 0;
        Objects.p62BodyPixels = 0; 
        Objects.p62BodyCount = 0; % Needed for invagination feature
        Objects.p62ShapeBySurface = 0; % Roundness feature
        Objects.p62BodycountByp62count = 0; % Invagination feature 
    end
    %% LAMP derived
      if sum(LAMPMask(:))> 0
        Objects.LAMPArea = sum(LAMPMask(:));
        Objects.LAMPAreaNorm = sum(LAMPMask(:))/sum(NucleiMask(:));
        voxelSizeX = 0.32;%Bin2
        voxelSizeY = 0.32; 
        Objects.MinLAMPVol = min(LAMPObjects.Area) * voxelSizeX * voxelSizeY;
        Objects.MaxLAMPVol = max(LAMPObjects.Area) * voxelSizeX * voxelSizeY;
        Objects.MeanLAMPVol = mean(LAMPObjects.Area) * voxelSizeX * voxelSizeY;
        Objects.StdLAMPVol = std(LAMPObjects.Area) * voxelSizeX * voxelSizeY;
        Objects.MedLAMPVol = median(LAMPObjects.Area) * voxelSizeX * voxelSizeY;
        Objects.MadLAMPVol = mad(LAMPObjects.Area, 1) * voxelSizeX * voxelSizeY;
        Objects.TotLAMPVolume = sum(LAMPMask(:)) * voxelSizeX * voxelSizeY;
        Objects.TotLAMPVolumeNorm = (sum(LAMPMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
        Objects.CountLAMP = size(LAMPObjects, 1);
        Objects.CountLAMPNorm = size(LAMPObjects, 1)/sum(NucleiMask(:));                
        Objects.MedMajorAxisLenghtLAMP = median(LAMPObjects.MajorAxisLength)* voxelSizeX;
        Objects.MedMinorAxisLenghtLAMP = median(LAMPObjects.MinorAxisLength)* voxelSizeX;
        Objects.MedLAMPIntensity = median(LAMPObjects.MeanIntensity);
        % Shape
        Conn6Strel = {};
        Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
        Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
        Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
        Conn6Strel = logical(cat(3, Conn6Strel{:}));
        LAMPErodedMask = imerode(LAMPMask, Conn6Strel);
        LAMPPerimMask = (LAMPMask - LAMPErodedMask) > 0;

        % Additional feature images and vectors
        LAMPBodyLabelIm = bwlabeln(LAMPErodedMask, 6);

        % Erosion derived
        Objects.LAMPPerimPixels = sum(LAMPPerimMask(:));
        Objects.LAMPBodyPixels = sum(LAMPErodedMask(:)); 
        Objects.LAMPBodyCount = max(LAMPBodyLabelIm(:)); % Needed for invagination feature
        Objects.LAMPShapeBySurface = Objects.LAMPBodyPixels / Objects.LAMPPerimPixels; % Roundness feature
        Objects.LAMPBodycountByLAMPcount = Objects.LAMPBodyCount / Objects.CountLAMP; % Invagination feature
      else
        Objects.LAMPArea = 0;
        Objects.LAMPAreaNorm = 0;
        Objects.MinLAMPVol = 0;
        Objects.MaxLAMPVol = 0;
        Objects.MeanLAMPVol = 0;
        Objects.StdLAMPVol = 0;
        Objects.MedLAMPVol = 0;
        Objects.MadLAMPVol = 0;
        Objects.TotLAMPVolume = 0;
        Objects.TotLAMPVolumeNorm = 0;
        Objects.CountLAMP = 0;
        Objects.CountLAMPNorm = 0; 
        Objects.MedMajorAxisLenghtLAMP = 0;
        Objects.MedMinorAxisLenghtLAMP = 0;
        Objects.MedLAMPIntensity =0;
        % Erosion derived
        Objects.LAMPPerimPixels = 0;
        Objects.LAMPBodyPixels = 0; 
        Objects.LAMPBodyCount = 0; % Needed for invagination feature
        Objects.LAMPShapeBySurface = 0; % Roundness feature
        Objects.LAMPBodycountByLAMPcount = 0; % Invagination feature 
      end
    
    %% ColoLp derived
    if sum(ColoLpMask(:))> 0
            Objects.ColoLpArea = sum(ColoLpMask(:));
            Objects.ColoLpAreaNorm = sum(ColoLpMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLpVol = min(ColoLpObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLpVol = max(ColoLpObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLpVol = mean(ColoLpObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLpVol = std(ColoLpObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLpVol = median(ColoLpObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLpVol = mad(ColoLpObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLpVolume = sum(ColoLpMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLpVolumeNorm = (sum(ColoLpMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLp = size(ColoLpObjects, 1);
            Objects.CountColoLpNorm = size(ColoLpObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLp = median(ColoLpObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLp = median(ColoLpObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLpErodedMask = imerode(ColoLpMask, Conn6Strel);
            ColoLpPerimMask = (ColoLpMask - ColoLpErodedMask) > 0;

            % Additional feature images and vectors
            ColoLpBodyLabelIm = bwlabeln(ColoLpErodedMask, 6);

            % Erosion derived
            Objects.ColoLpPerimPixels = sum(ColoLpPerimMask(:));
            Objects.ColoLpBodyPixels = sum(ColoLpErodedMask(:)); 
            Objects.ColoLpBodyCount = max(ColoLpBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLpShapeBySurface = Objects.ColoLpBodyPixels / Objects.ColoLpPerimPixels; % Roundness feature
            Objects.ColoLpBodycountByColoLpcount = Objects.ColoLpBodyCount / Objects.CountColoLp; % Invagination feature 
            % Previews
            ColoLpMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch3, [], 3)), zeros([size(ch1,1),size(ch1,2)], 'uint16'));
            ColoLpMaskPreview_perimeter = imoverlay(ColoLpMaskPreview, bwperim(max(ColoLpMask, [], 3)), [1 1 1]);
            %imtool(ColoLpMaskPreview)
            SavePathColoLpperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLp_peri.png'];
            imwrite(ColoLpMaskPreview_perimeter, SavePathColoLpperiPreview)
            else
            Objects.ColoLpArea = 0;
            Objects.ColoLpAreaNorm = 0;
            Objects.MinColoLpVol = 0;
            Objects.MaxColoLpVol = 0;
            Objects.MeanColoLpVol = 0;
            Objects.StdColoLpVol = 0;
            Objects.MedColoLpVol = 0;
            Objects.MadColoLpVol = 0;
            Objects.TotColoLpVolume = 0;
            Objects.TotColoLpVolumeNorm = 0;
            Objects.CountColoLp = 0;
            Objects.CountColoLpNorm = 0; 
            Objects.MedMajorAxisLenghtColoLp = 0;
            Objects.MedMinorAxisLenghtColoLp = 0;
    
			% Erosion derived
            Objects.ColoLpPerimPixels = 0;
            Objects.ColoLpBodyPixels = 0; 
            Objects.ColoLpBodyCount = 0; % Needed for invagination feature
            Objects.ColoLpShapeBySurface = 0; % Roundness feature
            Objects.ColoLpBodycountByColoLpcount = 0; % Invagination feature 
    end

%% ColoLpTS derived
    if sum(ColoLpTSMask(:))> 0
            Objects.ColoLpTSArea = sum(ColoLpTSMask(:));
            Objects.ColoLpTSAreaNorm = sum(ColoLpTSMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLpTSVol = min(ColoLpTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLpTSVol = max(ColoLpTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLpTSVol = mean(ColoLpTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLpTSVol = std(ColoLpTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLpTSVol = median(ColoLpTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLpTSVol = mad(ColoLpTSObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLpTSVolume = sum(ColoLpTSMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLpTSVolumeNorm = (sum(ColoLpTSMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLpTS = size(ColoLpTSObjects, 1);
            Objects.CountColoLpTSNorm = size(ColoLpTSObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLpTS = median(ColoLpTSObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLpTS = median(ColoLpTSObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLpTSErodedMask = imerode(ColoLpTSMask, Conn6Strel);
            ColoLpTSPerimMask = (ColoLpTSMask - ColoLpTSErodedMask) > 0;

            % Additional feature images and vectors
            ColoLpTSBodyLabelIm = bwlabeln(ColoLpTSErodedMask, 6);

            % Erosion derived
            Objects.ColoLpTSPerimPixels = sum(ColoLpTSPerimMask(:));
            Objects.ColoLpTSBodyPixels = sum(ColoLpTSErodedMask(:)); 
            Objects.ColoLpTSBodyCount = max(ColoLpTSBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLpTSShapeBySurface = Objects.ColoLpTSBodyPixels / Objects.ColoLpTSPerimPixels; % Roundness feature
            Objects.ColoLpTSBodycountByColoLpTScount = Objects.ColoLpTSBodyCount / Objects.CountColoLpTS; % Invagination feature 
            % Previews
            ColoLpTSMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch3, [], 3)), imadjust(max(ch4, [], 3)));
            ColoLpTSMaskPreview_perimeter = imoverlay(ColoLpTSMaskPreview, bwperim(max(ColoLpTSMask, [], 3)), [1 1 1]);
            %imtool(ColoLpTSMaskPreview)
            SavePathColoLpTSperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLpTS_peri.png'];
            imwrite(ColoLpTSMaskPreview_perimeter, SavePathColoLpTSperiPreview)
            else
            Objects.ColoLpTSArea = 0;
            Objects.ColoLpTSAreaNorm = 0;
            Objects.MinColoLpTSVol = 0;
            Objects.MaxColoLpTSVol = 0;
            Objects.MeanColoLpTSVol = 0;
            Objects.StdColoLpTSVol = 0;
            Objects.MedColoLpTSVol = 0;
            Objects.MadColoLpTSVol = 0;
            Objects.TotColoLpTSVolume = 0;
            Objects.TotColoLpTSVolumeNorm = 0;
            Objects.CountColoLpTS = 0;
            Objects.CountColoLpTSNorm = 0;                
			Objects.MedMajorAxisLenghtColoLpTS = 0;
            Objects.MedMinorAxisLenghtColoLpTS = 0;
    
            % Erosion derived
            Objects.ColoLpTSPerimPixels = 0;
            Objects.ColoLpTSBodyPixels = 0; 
            Objects.ColoLpTSBodyCount = 0; % Needed for invagination feature
            Objects.ColoLpTSShapeBySurface = 0; % Roundness feature
            Objects.ColoLpTSBodycountByColoLpTScount = 0; % Invagination feature 
    end

%% ColoLpTN derived
    if sum(ColoLpTNMask(:))> 0
            Objects.ColoLpTNArea = sum(ColoLpTNMask(:));
            Objects.ColoLpTNAreaNorm = sum(ColoLpTNMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLpTNVol = min(ColoLpTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLpTNVol = max(ColoLpTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLpTNVol = mean(ColoLpTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLpTNVol = std(ColoLpTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLpTNVol = median(ColoLpTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLpTNVol = mad(ColoLpTNObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLpTNVolume = sum(ColoLpTNMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLpTNVolumeNorm = (sum(ColoLpTNMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLpTN = size(ColoLpTNObjects, 1);
            Objects.CountColoLpTNNorm = size(ColoLpTNObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLpTN = median(ColoLpTNObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLpTN = median(ColoLpTNObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLpTNErodedMask = imerode(ColoLpTNMask, Conn6Strel);
            ColoLpTNPerimMask = (ColoLpTNMask - ColoLpTNErodedMask) > 0;

            % Additional feature images and vectors
            ColoLpTNBodyLabelIm = bwlabeln(ColoLpTNErodedMask, 6);

            % Erosion derived
            Objects.ColoLpTNPerimPixels = sum(ColoLpTNPerimMask(:));
            Objects.ColoLpTNBodyPixels = sum(ColoLpTNErodedMask(:)); 
            Objects.ColoLpTNBodyCount = max(ColoLpTNBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLpTNShapeBySurface = Objects.ColoLpTNBodyPixels / Objects.ColoLpTNPerimPixels; % Roundness feature
            Objects.ColoLpTNBodycountByColoLpTNcount = Objects.ColoLpTNBodyCount / Objects.CountColoLpTN; % Invagination feature 
            % Previews
            ColoLpTNMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch3, [], 3)), imadjust(max(ch4, [], 3)));
            ColoLpTNMaskPreview_perimeter = imoverlay(ColoLpTNMaskPreview, bwperim(max(ColoLpTNMask, [], 3)), [1 1 1]);
            %imtool(ColoLpTNMaskPreview)
            SavePathColoLpTNperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLpTN_peri.png'];
            imwrite(ColoLpTNMaskPreview_perimeter, SavePathColoLpTNperiPreview)
            else
            Objects.ColoLpTNArea = 0;
            Objects.ColoLpTNAreaNorm = 0;
            Objects.MinColoLpTNVol = 0;
            Objects.MaxColoLpTNVol = 0;
            Objects.MeanColoLpTNVol = 0;
            Objects.StdColoLpTNVol = 0;
            Objects.MedColoLpTNVol = 0;
            Objects.MadColoLpTNVol = 0;
            Objects.TotColoLpTNVolume = 0;
            Objects.TotColoLpTNVolumeNorm = 0;
            Objects.CountColoLpTN = 0;
            Objects.CountColoLpTNNorm = 0;                
			Objects.MedMajorAxisLenghtColoLpTN = 0;
            Objects.MedMinorAxisLenghtColoLpTN = 0;
    
            % Erosion derived
            Objects.ColoLpTNPerimPixels = 0;
            Objects.ColoLpTNBodyPixels = 0; 
            Objects.ColoLpTNBodyCount = 0; % Needed for invagination feature
            Objects.ColoLpTNShapeBySurface = 0; % Roundness feature
            Objects.ColoLpTNBodycountByColoLpTNcount = 0; % Invagination feature 
    end
%% ColopTS derived
    if sum(ColopTSMask(:))> 0
            Objects.ColopTSArea = sum(ColopTSMask(:));
            Objects.ColopTSAreaNorm = sum(ColopTSMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColopTSVol = min(ColopTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColopTSVol = max(ColopTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColopTSVol = mean(ColopTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColopTSVol = std(ColopTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColopTSVol = median(ColopTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColopTSVol = mad(ColopTSObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColopTSVolume = sum(ColopTSMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColopTSVolumeNorm = (sum(ColopTSMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColopTS = size(ColopTSObjects, 1);
            Objects.CountColopTSNorm = size(ColopTSObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColopTS = median(ColopTSObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColopTS = median(ColopTSObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColopTSErodedMask = imerode(ColopTSMask, Conn6Strel);
            ColopTSPerimMask = (ColopTSMask - ColopTSErodedMask) > 0;

            % Additional feature images and vectors
            ColopTSBodyLabelIm = bwlabeln(ColopTSErodedMask, 6);

            % Erosion derived
            Objects.ColopTSPerimPixels = sum(ColopTSPerimMask(:));
            Objects.ColopTSBodyPixels = sum(ColopTSErodedMask(:)); 
            Objects.ColopTSBodyCount = max(ColopTSBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColopTSShapeBySurface = Objects.ColopTSBodyPixels / Objects.ColopTSPerimPixels; % Roundness feature
            Objects.ColopTSBodycountByColopTScount = Objects.ColopTSBodyCount / Objects.CountColopTS; % Invagination feature 
            % Previews
            ColopTSMaskPreview = cat(3, imadjust(max(ch3, [], 3)), imadjust(max(ch4, [], 3)), imadjust(max(ch1, [], 3)));
            ColopTSMaskPreview_perimeter = imoverlay(ColopTSMaskPreview, bwperim(max(ColopTSMask, [], 3)), [1 1 1]);
            %imtool(ColopTSMaskPreview)
            SavePathColopTSperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColopTS_peri.png'];
            imwrite(ColopTSMaskPreview_perimeter, SavePathColopTSperiPreview)
            else
            Objects.ColopTSArea = 0;
            Objects.ColopTSAreaNorm = 0;
            Objects.MinColopTSVol = 0;
            Objects.MaxColopTSVol = 0;
            Objects.MeanColopTSVol = 0;
            Objects.StdColopTSVol = 0;
            Objects.MedColopTSVol = 0;
            Objects.MadColopTSVol = 0;
            Objects.TotColopTSVolume = 0;
            Objects.TotColopTSVolumeNorm = 0;
            Objects.CountColopTS = 0;
            Objects.CountColopTSNorm = 0;                
			Objects.MedMajorAxisLenghtColopTS = 0;
            Objects.MedMinorAxisLenghtColopTS = 0;
    
            % Erosion derived
            Objects.ColopTSPerimPixels = 0;
            Objects.ColopTSBodyPixels = 0; 
            Objects.ColopTSBodyCount = 0; % Needed for invagination feature
            Objects.ColopTSShapeBySurface = 0; % Roundness feature
            Objects.ColopTSBodycountByColopTScount = 0; % Invagination feature 
    end
%% ColopTN derived
    if sum(ColopTNMask(:))> 0
            Objects.ColopTNArea = sum(ColopTNMask(:));
            Objects.ColopTNAreaNorm = sum(ColopTNMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColopTNVol = min(ColopTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColopTNVol = max(ColopTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColopTNVol = mean(ColopTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColopTNVol = std(ColopTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColopTNVol = median(ColopTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColopTNVol = mad(ColopTNObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColopTNVolume = sum(ColopTNMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColopTNVolumeNorm = (sum(ColopTNMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColopTN = size(ColopTNObjects, 1);
            Objects.CountColopTNNorm = size(ColopTNObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColopTN = median(ColopTNObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColopTN = median(ColopTNObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColopTNErodedMask = imerode(ColopTNMask, Conn6Strel);
            ColopTNPerimMask = (ColopTNMask - ColopTNErodedMask) > 0;

            % Additional feature images and vectors
            ColopTNBodyLabelIm = bwlabeln(ColopTNErodedMask, 6);

            % Erosion derived
            Objects.ColopTNPerimPixels = sum(ColopTNPerimMask(:));
            Objects.ColopTNBodyPixels = sum(ColopTNErodedMask(:)); 
            Objects.ColopTNBodyCount = max(ColopTNBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColopTNShapeBySurface = Objects.ColopTNBodyPixels / Objects.ColopTNPerimPixels; % Roundness feature
            Objects.ColopTNBodycountByColopTNcount = Objects.ColopTNBodyCount / Objects.CountColopTN; % Invagination feature 
            % Previews
            ColopTNMaskPreview = cat(3, imadjust(max(ch3, [], 3)), imadjust(max(ch4, [], 3)), imadjust(max(ch1, [], 3)));
            ColopTNMaskPreview_perimeter = imoverlay(ColopTNMaskPreview, bwperim(max(ColopTNMask, [], 3)), [1 1 1]);
            %imtool(ColopTNMaskPreview)
            SavePathColopTNperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColopTN_peri.png'];
            imwrite(ColopTNMaskPreview_perimeter, SavePathColopTNperiPreview)
            else
            Objects.ColopTNArea = 0;
            Objects.ColopTNAreaNorm = 0;
            Objects.MinColopTNVol = 0;
            Objects.MaxColopTNVol = 0;
            Objects.MeanColopTNVol = 0;
            Objects.StdColopTNVol = 0;
            Objects.MedColopTNVol = 0;
            Objects.MadColopTNVol = 0;
            Objects.TotColopTNVolume = 0;
            Objects.TotColopTNVolumeNorm = 0;
            Objects.CountColopTN = 0;
            Objects.CountColopTNNorm = 0;                
			Objects.MedMajorAxisLenghtColopTN = 0;
            Objects.MedMinorAxisLenghtColopTN = 0;
    
            % Erosion derived
            Objects.ColopTNPerimPixels = 0;
            Objects.ColopTNBodyPixels = 0; 
            Objects.ColopTNBodyCount = 0; % Needed for invagination feature
            Objects.ColopTNShapeBySurface = 0; % Roundness feature
            Objects.ColopTNBodycountByColopTNcount = 0; % Invagination feature 
    end
%% ColoLATS derived
    if sum(ColoLATSMask(:))> 0
            Objects.ColoLATSArea = sum(ColoLATSMask(:));
            Objects.ColoLATSAreaNorm = sum(ColoLATSMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLATSVol = min(ColoLATSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLATSVol = max(ColoLATSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLATSVol = mean(ColoLATSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLATSVol = std(ColoLATSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLATSVol = median(ColoLATSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLATSVol = mad(ColoLATSObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLATSVolume = sum(ColoLATSMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLATSVolumeNorm = (sum(ColoLATSMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLATS = size(ColoLATSObjects, 1);
            Objects.CountColoLATSNorm = size(ColoLATSObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLATS = median(ColoLATSObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLATS = median(ColoLATSObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLATSErodedMask = imerode(ColoLATSMask, Conn6Strel);
            ColoLATSPerimMask = (ColoLATSMask - ColoLATSErodedMask) > 0;

            % Additional feature images and vectors
            ColoLATSBodyLabelIm = bwlabeln(ColoLATSErodedMask, 6);

            % Erosion derived
            Objects.ColoLATSPerimPixels = sum(ColoLATSPerimMask(:));
            Objects.ColoLATSBodyPixels = sum(ColoLATSErodedMask(:)); 
            Objects.ColoLATSBodyCount = max(ColoLATSBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLATSShapeBySurface = Objects.ColoLATSBodyPixels / Objects.ColoLATSPerimPixels; % Roundness feature
            Objects.ColoLATSBodycountByColoLATScount = Objects.ColoLATSBodyCount / Objects.CountColoLATS; % Invagination feature 
            % Previews
            ColoLATSMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch4, [], 3)), imadjust(max(ch1, [], 3)));
            ColoLATSMaskPreview_perimeter = imoverlay(ColoLATSMaskPreview, bwperim(max(ColoLATSMask, [], 3)), [1 1 1]);
            %imtool(ColoLATSMaskPreview)
            SavePathColoLATSperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLATS_peri.png'];
            imwrite(ColoLATSMaskPreview_perimeter, SavePathColoLATSperiPreview)
            else
            Objects.ColoLATSArea = 0;
            Objects.ColoLATSAreaNorm = 0;
            Objects.MinColoLATSVol = 0;
            Objects.MaxColoLATSVol = 0;
            Objects.MeanColoLATSVol = 0;
            Objects.StdColoLATSVol = 0;
            Objects.MedColoLATSVol = 0;
            Objects.MadColoLATSVol = 0;
            Objects.TotColoLATSVolume = 0;
            Objects.TotColoLATSVolumeNorm = 0;
            Objects.CountColoLATS = 0;
            Objects.CountColoLATSNorm = 0;                
			Objects.MedMajorAxisLenghtColoLATS = 0;
            Objects.MedMinorAxisLenghtColoLATS = 0;
    
            % Erosion derived
            Objects.ColoLATSPerimPixels = 0;
            Objects.ColoLATSBodyPixels = 0; 
            Objects.ColoLATSBodyCount = 0; % Needed for invagination feature
            Objects.ColoLATSShapeBySurface = 0; % Roundness feature
            Objects.ColoLATSBodycountByColoLATScount = 0; % Invagination feature 
    end
%% ColoLATN derived
    if sum(ColoLATNMask(:))> 0
            Objects.ColoLATNArea = sum(ColoLATNMask(:));
            Objects.ColoLATNAreaNorm = sum(ColoLATNMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLATNVol = min(ColoLATNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLATNVol = max(ColoLATNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLATNVol = mean(ColoLATNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLATNVol = std(ColoLATNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLATNVol = median(ColoLATNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLATNVol = mad(ColoLATNObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLATNVolume = sum(ColoLATNMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLATNVolumeNorm = (sum(ColoLATNMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLATN = size(ColoLATNObjects, 1);
            Objects.CountColoLATNNorm = size(ColoLATNObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLATN = median(ColoLATNObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLATN = median(ColoLATNObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLATNErodedMask = imerode(ColoLATNMask, Conn6Strel);
            ColoLATNPerimMask = (ColoLATNMask - ColoLATNErodedMask) > 0;

            % Additional feature images and vectors
            ColoLATNBodyLabelIm = bwlabeln(ColoLATNErodedMask, 6);

            % Erosion derived
            Objects.ColoLATNPerimPixels = sum(ColoLATNPerimMask(:));
            Objects.ColoLATNBodyPixels = sum(ColoLATNErodedMask(:)); 
            Objects.ColoLATNBodyCount = max(ColoLATNBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLATNShapeBySurface = Objects.ColoLATNBodyPixels / Objects.ColoLATNPerimPixels; % Roundness feature
            Objects.ColoLATNBodycountByColoLATNcount = Objects.ColoLATNBodyCount / Objects.CountColoLATN; % Invagination feature 
            % Previews
            ColoLATNMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch4, [], 3)), imadjust(max(ch1, [], 3)));
            ColoLATNMaskPreview_perimeter = imoverlay(ColoLATNMaskPreview, bwperim(max(ColoLATNMask, [], 3)), [1 1 1]);
            %imtool(ColoLATNMaskPreview)
            SavePathColoLATNperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLATN_peri.png'];
            imwrite(ColoLATNMaskPreview_perimeter, SavePathColoLATNperiPreview)
            else
            Objects.ColoLATNArea = 0;
            Objects.ColoLATNAreaNorm = 0;
            Objects.MinColoLATNVol = 0;
            Objects.MaxColoLATNVol = 0;
            Objects.MeanColoLATNVol = 0;
            Objects.StdColoLATNVol = 0;
            Objects.MedColoLATNVol = 0;
            Objects.MadColoLATNVol = 0;
            Objects.TotColoLATNVolume = 0;
            Objects.TotColoLATNVolumeNorm = 0;
            Objects.CountColoLATN = 0;
            Objects.CountColoLATNNorm = 0;                
			Objects.MedMajorAxisLenghtColoLATN = 0;
            Objects.MedMinorAxisLenghtColoLATN = 0;
    
            % Erosion derived
            Objects.ColoLATNPerimPixels = 0;
            Objects.ColoLATNBodyPixels = 0; 
            Objects.ColoLATNBodyCount = 0; % Needed for invagination feature
            Objects.ColoLATNShapeBySurface = 0; % Roundness feature
            Objects.ColoLATNBodycountByColoLATNcount = 0; % Invagination feature 
    end	
end

