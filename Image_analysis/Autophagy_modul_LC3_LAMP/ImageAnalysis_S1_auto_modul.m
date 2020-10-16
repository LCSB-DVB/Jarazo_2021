function [Objects] = ImageAnalysis_S1_auto_modul(ch1,ch2,ch3,ch4,InfoTableThis, PreviewPath)
  
%% Segment nuclei
    NucleiBlurred = imfilter(ch1, fspecial('gaussian',5, 4)); % it(NucleiBlurred)
    NucleiMask = NucleiBlurred > 120; % it(NucleiMask)
        if sum(NucleiMask(:))== 0
            Objects = {}; 
            return
        end
    NucleiDeadMask = NucleiBlurred > 1500; % it(NucleiDeadMask)
    %% Segment LC3somes
    LC3DoG = imfilter(ch2, fspecial('gaussian', 5, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 5, 2), 'symmetric'); % it(LC3DoG)
    LC3MaskLocal = LC3DoG > 10; % it(LC3MaskLocal)
    LC3Mask = LC3MaskLocal; % & LC3GlobalMask;
    LC3Mask = bwareaopen (LC3Mask, 2); % it(LC3Mask)
   
    %% Segment LAMPsomes
    LAMPDoG = imfilter(ch3, fspecial('gaussian', 5, 1), 'symmetric') - imfilter(ch3, fspecial('gaussian', 5, 2), 'symmetric'); % it(LAMPDoG)
    LAMPMaskLocal = LAMPDoG > 6; % it(LAMPMaskLocal)
    LAMPMask = LAMPMaskLocal; % & LAMPGlobalMask;
    LAMPMask = bwareaopen (LAMPMask, 2); % it(LAMPMask)
        
    %% Segment TH Soma
    THBlurred = imfilter(ch4, fspecial('gaussian',20, 4)); % it(THBlurred)
    THSomaMask = THBlurred > 100;  % it(THSomaMask)
    THSomaMask = bwareaopen(THSomaMask, 600);
	
	%% Segment TH Neurite
	THNeuriteDoG = imfilter(ch4, fspecial('gaussian', 7, 1), 'symmetric') - imfilter(ch4, fspecial('gaussian', 7, 2), 'symmetric');  %it(THNeuriteDoG)
    THNeuriteMask = THNeuriteDoG > 5;  % it(THNeuriteMask)
    THNeuriteMask = bwareaopen(THNeuriteMask, 7);
    THNeuriteMask = THNeuriteMask & ~THSomaMask;
	
	%% LC3-LAMP Colocalization Mask
    ColoLLMask = LC3Mask & LAMPMask; %it(ColoLLMask) 
	    
	%% LC3-LAMP-THSoma Colocalization Mask
    ColoLLTSMask = LC3Mask & LAMPMask & THSomaMask; %it(ColoLLTSMask) 
    
	%% LC3-LAMP-THNeurite Colocalization Mask
    ColoLLTNMask = LC3Mask & LAMPMask & THNeuriteMask; %it(ColoLLTNMask) 
	
	%% LC3-THSoma Colocalization Mask
    ColoLCTSMask = LC3Mask & THSomaMask; %it(ColoLCTSMask) 

	%% LC3-THSoma Colocalization Mask
    ColoLCTNMask = LC3Mask & THNeuriteMask; %it(ColoLCTNMask) 

	%% LAMP-THSoma Colocalization Mask
    ColoLATSMask = LAMPMask & THSomaMask; %it(ColoLATSMask) 
	
	%% LAMP-THNeurite Colocalization Mask
    ColoLATNMask = LAMPMask & THNeuriteMask; %it(ColoLATNMask) 
	
	%% Morphometrics
    LC3LM = bwlabeln(LC3Mask);
    LC3Objects = regionprops('table', LC3LM, ch2,{'Area','MajorAxisLength','MinorAxisLength','MeanIntensity'});
    LAMPLM = bwlabeln(LAMPMask);
    LAMPObjects = regionprops('table', LAMPLM, ch3,{'Area','MajorAxisLength','MinorAxisLength','MeanIntensity'});
    ColoLLLM = bwlabeln(ColoLLMask);
    ColoLLObjects = regionprops('table', ColoLLLM,{'Area','MajorAxisLength','MinorAxisLength'});
    ColoLLTSLM = bwlabeln(ColoLLTSMask);
    ColoLLTSObjects = regionprops('table', ColoLLTSLM, {'Area','MajorAxisLength','MinorAxisLength'});
	ColoLLTNLM = bwlabeln(ColoLLTNMask);
    ColoLLTNObjects = regionprops('table', ColoLLTNLM, {'Area','MajorAxisLength','MinorAxisLength'});
    ColoLCTSLM = bwlabeln(ColoLCTSMask);
    ColoLCTSObjects = regionprops('table', ColoLCTSLM, {'Area','MajorAxisLength','MinorAxisLength'});
	ColoLCTNLM = bwlabeln(ColoLCTNMask);
    ColoLCTNObjects = regionprops('table', ColoLCTNLM, {'Area','MajorAxisLength','MinorAxisLength'});
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
    LC3MaskPreview = imoverlay(imadjust(max(ch2, [], 3)), bwperim(max(LC3Mask, [], 3)), [1 0 0]);
    LC3MaskPreview = imoverlay(LC3MaskPreview, BarMask, [1 1 1]);
    %imtool(LC3MaskPreview )
    LC3RawPreview = imadjust(max(ch2, [], 3));
    LC3RawPreview = imoverlay(LC3RawPreview, BarMask, [1 1 1]);
    %imtool(LC3RawPreview)
    
    LAMPMaskPreview = imoverlay(imadjust(max(ch3, [], 3)), bwperim(max(LAMPMask, [], 3)), [1 0 0]);
    LAMPMaskPreview = imoverlay(LAMPMaskPreview, BarMask, [1 1 1]);
    %imtool(LAMPMaskPreview )
    LAMPRawPreview = imadjust(max(ch3, [], 3));
    LAMPRawPreview = imoverlay(LAMPRawPreview, BarMask, [1 1 1]);
    %imtool(LAMPRawPreview)
    
   
    SavePathLC3MaskPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_LC3_mask.png'];
    SavePathLC3rawPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_LC3_raw.png'];
    SavePathLAMPMaskPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_LAMP_mask.png'];
    SavePathLAMPrawPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_LAMP_raw.png'];
    SavePathRGBPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RGB.png'];
    SavePathRGB_1Preview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RGB_1.png'];
    
   
    imwrite(LC3MaskPreview, SavePathLC3MaskPreview)
    imwrite(LC3RawPreview,  SavePathLC3rawPreview)
    imwrite(LAMPMaskPreview, SavePathLAMPMaskPreview)
    imwrite(LAMPRawPreview,  SavePathLAMPrawPreview)
    imwrite(RGB, SavePathRGBPreview)
    imwrite(RGB_1, SavePathRGB_1Preview)
    
    %% LC3 derived
    if sum(LC3Mask(:))> 0
        Objects.LC3Area = sum(LC3Mask(:));
        Objects.LC3AreaNorm = sum(LC3Mask(:))/sum(NucleiMask(:));
        voxelSizeX = 0.32;%Bin2
        voxelSizeY = 0.32; 
        Objects.MinLC3Vol = min(LC3Objects.Area) * voxelSizeX * voxelSizeY; 
        Objects.MaxLC3Vol = max(LC3Objects.Area) * voxelSizeX * voxelSizeY;
        Objects.MeanLC3Vol = mean(LC3Objects.Area) * voxelSizeX * voxelSizeY ;
        Objects.StdLC3Vol = std(LC3Objects.Area) * voxelSizeX * voxelSizeY;
        Objects.MedLC3Vol = median(LC3Objects.Area) * voxelSizeX * voxelSizeY;
        Objects.MadLC3Vol = mad(LC3Objects.Area, 1) * voxelSizeX * voxelSizeY;
        Objects.TotLC3Volume = sum(LC3Mask(:)) * voxelSizeX * voxelSizeY;
        Objects.TotLC3VolumeNorm = (sum(LC3Mask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
        Objects.CountLC3 = size(LC3Objects, 1);
        Objects.CountLC3Norm = size(LC3Objects, 1)/sum(NucleiMask(:));                
        Objects.MedMajorAxisLenghtLC3 = median(LC3Objects.MajorAxisLength)* voxelSizeX;
        Objects.MedMinorAxisLenghtLC3 = median(LC3Objects.MinorAxisLength)* voxelSizeX;
        Objects.MedLC3Intensity = median(LC3Objects.MeanIntensity);
        % Shape
        Conn6Strel = {};
        Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
        Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
        Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
        Conn6Strel = logical(cat(3, Conn6Strel{:}));
        LC3ErodedMask = imerode(LC3Mask, Conn6Strel);
        LC3PerimMask = (LC3Mask - LC3ErodedMask) > 0;

        % Additional feature images and vectors
        LC3BodyLabelIm = bwlabeln(LC3ErodedMask, 6);

        % Erosion derived
        Objects.LC3PerimPixels = sum(LC3PerimMask(:));
        Objects.LC3BodyPixels = sum(LC3ErodedMask(:)); 
        Objects.LC3BodyCount = max(LC3BodyLabelIm(:)); % Needed for invagination feature
        Objects.LC3ShapeBySurface = Objects.LC3BodyPixels / Objects.LC3PerimPixels; % Roundness feature
        Objects.LC3BodycountByLC3count = Objects.LC3BodyCount / Objects.CountLC3; % Invagination feature
        else
        Objects.LC3Area = 0;
        Objects.LC3AreaNorm = 0;
        Objects.MinLC3Vol = 0;
        Objects.MaxLC3Vol = 0;
        Objects.MeanLC3Vol = 0;
        Objects.StdLC3Vol = 0;
        Objects.MedLC3Vol = 0;
        Objects.MadLC3Vol = 0;
        Objects.TotLC3Volume = 0;
        Objects.TotLC3VolumeNorm = 0;
        Objects.CountLC3 = 0;
        Objects.CountLC3Norm = 0; 
        Objects.MedMajorAxisLenghtLC3 = 0;
        Objects.MedMinorAxisLenghtLC3 = 0;
        Objects.MedLC3Intensity =0;
        % Erosion derived
        Objects.LC3PerimPixels = 0;
        Objects.LC3BodyPixels = 0; 
        Objects.LC3BodyCount = 0; % Needed for invagination feature
        Objects.LC3ShapeBySurface = 0; % Roundness feature
        Objects.LC3BodycountByLC3count = 0; % Invagination feature 
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
    
    %% ColoLL derived
    if sum(ColoLLMask(:))> 0
            Objects.ColoLLArea = sum(ColoLLMask(:));
            Objects.ColoLLAreaNorm = sum(ColoLLMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLLVol = min(ColoLLObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLLVol = max(ColoLLObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLLVol = mean(ColoLLObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLLVol = std(ColoLLObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLLVol = median(ColoLLObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLLVol = mad(ColoLLObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLLVolume = sum(ColoLLMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLLVolumeNorm = (sum(ColoLLMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLL = size(ColoLLObjects, 1);
            Objects.CountColoLLNorm = size(ColoLLObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLL = median(ColoLLObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLL = median(ColoLLObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLLErodedMask = imerode(ColoLLMask, Conn6Strel);
            ColoLLPerimMask = (ColoLLMask - ColoLLErodedMask) > 0;

            % Additional feature images and vectors
            ColoLLBodyLabelIm = bwlabeln(ColoLLErodedMask, 6);

            % Erosion derived
            Objects.ColoLLPerimPixels = sum(ColoLLPerimMask(:));
            Objects.ColoLLBodyPixels = sum(ColoLLErodedMask(:)); 
            Objects.ColoLLBodyCount = max(ColoLLBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLLShapeBySurface = Objects.ColoLLBodyPixels / Objects.ColoLLPerimPixels; % Roundness feature
            Objects.ColoLLBodycountByColoLLcount = Objects.ColoLLBodyCount / Objects.CountColoLL; % Invagination feature 
            % Previews
            ColoLLMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch3, [], 3)), zeros([size(ch1,1),size(ch1,2)], 'uint16'));
            ColoLLMaskPreview_perimeter = imoverlay(ColoLLMaskPreview, bwperim(max(ColoLLMask, [], 3)), [1 1 1]);
            %imtool(ColoLLMaskPreview)
            SavePathColoLLperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLL_peri.png'];
            imwrite(ColoLLMaskPreview_perimeter, SavePathColoLLperiPreview)
            else
            Objects.ColoLLArea = 0;
            Objects.ColoLLAreaNorm = 0;
            Objects.MinColoLLVol = 0;
            Objects.MaxColoLLVol = 0;
            Objects.MeanColoLLVol = 0;
            Objects.StdColoLLVol = 0;
            Objects.MedColoLLVol = 0;
            Objects.MadColoLLVol = 0;
            Objects.TotColoLLVolume = 0;
            Objects.TotColoLLVolumeNorm = 0;
            Objects.CountColoLL = 0;
            Objects.CountColoLLNorm = 0; 
            Objects.MedMajorAxisLenghtColoLL = 0;
            Objects.MedMinorAxisLenghtColoLL = 0;
    
			% Erosion derived
            Objects.ColoLLPerimPixels = 0;
            Objects.ColoLLBodyPixels = 0; 
            Objects.ColoLLBodyCount = 0; % Needed for invagination feature
            Objects.ColoLLShapeBySurface = 0; % Roundness feature
            Objects.ColoLLBodycountByColoLLcount = 0; % Invagination feature 
    end

%% ColoLLTS derived
    if sum(ColoLLTSMask(:))> 0
            Objects.ColoLLTSArea = sum(ColoLLTSMask(:));
            Objects.ColoLLTSAreaNorm = sum(ColoLLTSMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLLTSVol = min(ColoLLTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLLTSVol = max(ColoLLTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLLTSVol = mean(ColoLLTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLLTSVol = std(ColoLLTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLLTSVol = median(ColoLLTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLLTSVol = mad(ColoLLTSObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLLTSVolume = sum(ColoLLTSMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLLTSVolumeNorm = (sum(ColoLLTSMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLLTS = size(ColoLLTSObjects, 1);
            Objects.CountColoLLTSNorm = size(ColoLLTSObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLLTS = median(ColoLLTSObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLLTS = median(ColoLLTSObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLLTSErodedMask = imerode(ColoLLTSMask, Conn6Strel);
            ColoLLTSPerimMask = (ColoLLTSMask - ColoLLTSErodedMask) > 0;

            % Additional feature images and vectors
            ColoLLTSBodyLabelIm = bwlabeln(ColoLLTSErodedMask, 6);

            % Erosion derived
            Objects.ColoLLTSPerimPixels = sum(ColoLLTSPerimMask(:));
            Objects.ColoLLTSBodyPixels = sum(ColoLLTSErodedMask(:)); 
            Objects.ColoLLTSBodyCount = max(ColoLLTSBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLLTSShapeBySurface = Objects.ColoLLTSBodyPixels / Objects.ColoLLTSPerimPixels; % Roundness feature
            Objects.ColoLLTSBodycountByColoLLTScount = Objects.ColoLLTSBodyCount / Objects.CountColoLLTS; % Invagination feature 
            % Previews
            ColoLLTSMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch3, [], 3)), imadjust(max(ch4, [], 3)));
            ColoLLTSMaskPreview_perimeter = imoverlay(ColoLLTSMaskPreview, bwperim(max(ColoLLTSMask, [], 3)), [1 1 1]);
            %imtool(ColoLLTSMaskPreview)
            SavePathColoLLTSperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLLTS_peri.png'];
            imwrite(ColoLLTSMaskPreview_perimeter, SavePathColoLLTSperiPreview)
            else
            Objects.ColoLLTSArea = 0;
            Objects.ColoLLTSAreaNorm = 0;
            Objects.MinColoLLTSVol = 0;
            Objects.MaxColoLLTSVol = 0;
            Objects.MeanColoLLTSVol = 0;
            Objects.StdColoLLTSVol = 0;
            Objects.MedColoLLTSVol = 0;
            Objects.MadColoLLTSVol = 0;
            Objects.TotColoLLTSVolume = 0;
            Objects.TotColoLLTSVolumeNorm = 0;
            Objects.CountColoLLTS = 0;
            Objects.CountColoLLTSNorm = 0;                
			Objects.MedMajorAxisLenghtColoLLTS = 0;
            Objects.MedMinorAxisLenghtColoLLTS = 0;
    
            % Erosion derived
            Objects.ColoLLTSPerimPixels = 0;
            Objects.ColoLLTSBodyPixels = 0; 
            Objects.ColoLLTSBodyCount = 0; % Needed for invagination feature
            Objects.ColoLLTSShapeBySurface = 0; % Roundness feature
            Objects.ColoLLTSBodycountByColoLLTScount = 0; % Invagination feature 
    end

%% ColoLLTN derived
    if sum(ColoLLTNMask(:))> 0
            Objects.ColoLLTNArea = sum(ColoLLTNMask(:));
            Objects.ColoLLTNAreaNorm = sum(ColoLLTNMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLLTNVol = min(ColoLLTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLLTNVol = max(ColoLLTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLLTNVol = mean(ColoLLTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLLTNVol = std(ColoLLTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLLTNVol = median(ColoLLTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLLTNVol = mad(ColoLLTNObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLLTNVolume = sum(ColoLLTNMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLLTNVolumeNorm = (sum(ColoLLTNMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLLTN = size(ColoLLTNObjects, 1);
            Objects.CountColoLLTNNorm = size(ColoLLTNObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLLTN = median(ColoLLTNObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLLTN = median(ColoLLTNObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLLTNErodedMask = imerode(ColoLLTNMask, Conn6Strel);
            ColoLLTNPerimMask = (ColoLLTNMask - ColoLLTNErodedMask) > 0;

            % Additional feature images and vectors
            ColoLLTNBodyLabelIm = bwlabeln(ColoLLTNErodedMask, 6);

            % Erosion derived
            Objects.ColoLLTNPerimPixels = sum(ColoLLTNPerimMask(:));
            Objects.ColoLLTNBodyPixels = sum(ColoLLTNErodedMask(:)); 
            Objects.ColoLLTNBodyCount = max(ColoLLTNBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLLTNShapeBySurface = Objects.ColoLLTNBodyPixels / Objects.ColoLLTNPerimPixels; % Roundness feature
            Objects.ColoLLTNBodycountByColoLLTNcount = Objects.ColoLLTNBodyCount / Objects.CountColoLLTN; % Invagination feature 
            % Previews
            ColoLLTNMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch3, [], 3)), imadjust(max(ch4, [], 3)));
            ColoLLTNMaskPreview_perimeter = imoverlay(ColoLLTNMaskPreview, bwperim(max(ColoLLTNMask, [], 3)), [1 1 1]);
            %imtool(ColoLLTNMaskPreview)
            SavePathColoLLTNperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLLTN_peri.png'];
            imwrite(ColoLLTNMaskPreview_perimeter, SavePathColoLLTNperiPreview)
            else
            Objects.ColoLLTNArea = 0;
            Objects.ColoLLTNAreaNorm = 0;
            Objects.MinColoLLTNVol = 0;
            Objects.MaxColoLLTNVol = 0;
            Objects.MeanColoLLTNVol = 0;
            Objects.StdColoLLTNVol = 0;
            Objects.MedColoLLTNVol = 0;
            Objects.MadColoLLTNVol = 0;
            Objects.TotColoLLTNVolume = 0;
            Objects.TotColoLLTNVolumeNorm = 0;
            Objects.CountColoLLTN = 0;
            Objects.CountColoLLTNNorm = 0;                
			Objects.MedMajorAxisLenghtColoLLTN = 0;
            Objects.MedMinorAxisLenghtColoLLTN = 0;
    
            % Erosion derived
            Objects.ColoLLTNPerimPixels = 0;
            Objects.ColoLLTNBodyPixels = 0; 
            Objects.ColoLLTNBodyCount = 0; % Needed for invagination feature
            Objects.ColoLLTNShapeBySurface = 0; % Roundness feature
            Objects.ColoLLTNBodycountByColoLLTNcount = 0; % Invagination feature 
    end
%% ColoLCTS derived
    if sum(ColoLCTSMask(:))> 0
            Objects.ColoLCTSArea = sum(ColoLCTSMask(:));
            Objects.ColoLCTSAreaNorm = sum(ColoLCTSMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLCTSVol = min(ColoLCTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLCTSVol = max(ColoLCTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLCTSVol = mean(ColoLCTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLCTSVol = std(ColoLCTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLCTSVol = median(ColoLCTSObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLCTSVol = mad(ColoLCTSObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLCTSVolume = sum(ColoLCTSMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLCTSVolumeNorm = (sum(ColoLCTSMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLCTS = size(ColoLCTSObjects, 1);
            Objects.CountColoLCTSNorm = size(ColoLCTSObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLCTS = median(ColoLCTSObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLCTS = median(ColoLCTSObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLCTSErodedMask = imerode(ColoLCTSMask, Conn6Strel);
            ColoLCTSPerimMask = (ColoLCTSMask - ColoLCTSErodedMask) > 0;

            % Additional feature images and vectors
            ColoLCTSBodyLabelIm = bwlabeln(ColoLCTSErodedMask, 6);

            % Erosion derived
            Objects.ColoLCTSPerimPixels = sum(ColoLCTSPerimMask(:));
            Objects.ColoLCTSBodyPixels = sum(ColoLCTSErodedMask(:)); 
            Objects.ColoLCTSBodyCount = max(ColoLCTSBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLCTSShapeBySurface = Objects.ColoLCTSBodyPixels / Objects.ColoLCTSPerimPixels; % Roundness feature
            Objects.ColoLCTSBodycountByColoLCTScount = Objects.ColoLCTSBodyCount / Objects.CountColoLCTS; % Invagination feature 
            % Previews
            ColoLCTSMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch4, [], 3)), imadjust(max(ch1, [], 3)));
            ColoLCTSMaskPreview_perimeter = imoverlay(ColoLCTSMaskPreview, bwperim(max(ColoLCTSMask, [], 3)), [1 1 1]);
            %imtool(ColoLCTSMaskPreview)
            SavePathColoLCTSperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLCTS_peri.png'];
            imwrite(ColoLCTSMaskPreview_perimeter, SavePathColoLCTSperiPreview)
            else
            Objects.ColoLCTSArea = 0;
            Objects.ColoLCTSAreaNorm = 0;
            Objects.MinColoLCTSVol = 0;
            Objects.MaxColoLCTSVol = 0;
            Objects.MeanColoLCTSVol = 0;
            Objects.StdColoLCTSVol = 0;
            Objects.MedColoLCTSVol = 0;
            Objects.MadColoLCTSVol = 0;
            Objects.TotColoLCTSVolume = 0;
            Objects.TotColoLCTSVolumeNorm = 0;
            Objects.CountColoLCTS = 0;
            Objects.CountColoLCTSNorm = 0;                
			Objects.MedMajorAxisLenghtColoLCTS = 0;
            Objects.MedMinorAxisLenghtColoLCTS = 0;
    
            % Erosion derived
            Objects.ColoLCTSPerimPixels = 0;
            Objects.ColoLCTSBodyPixels = 0; 
            Objects.ColoLCTSBodyCount = 0; % Needed for invagination feature
            Objects.ColoLCTSShapeBySurface = 0; % Roundness feature
            Objects.ColoLCTSBodycountByColoLCTScount = 0; % Invagination feature 
    end
%% ColoLCTN derived
    if sum(ColoLCTNMask(:))> 0
            Objects.ColoLCTNArea = sum(ColoLCTNMask(:));
            Objects.ColoLCTNAreaNorm = sum(ColoLCTNMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.32;%Bin2
            voxelSizeY = 0.32; 
            Objects.MinColoLCTNVol = min(ColoLCTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MaxColoLCTNVol = max(ColoLCTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MeanColoLCTNVol = mean(ColoLCTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.StdColoLCTNVol = std(ColoLCTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MedColoLCTNVol = median(ColoLCTNObjects.Area) * voxelSizeX * voxelSizeY;
            Objects.MadColoLCTNVol = mad(ColoLCTNObjects.Area, 1) * voxelSizeX * voxelSizeY;
            Objects.TotColoLCTNVolume = sum(ColoLCTNMask(:)) * voxelSizeX * voxelSizeY;
            Objects.TotColoLCTNVolumeNorm = (sum(ColoLCTNMask(:)) * voxelSizeX * voxelSizeY)/sum(NucleiMask(:));
            Objects.CountColoLCTN = size(ColoLCTNObjects, 1);
            Objects.CountColoLCTNNorm = size(ColoLCTNObjects, 1)/sum(NucleiMask(:));                
            Objects.MedMajorAxisLenghtColoLCTN = median(ColoLCTNObjects.MajorAxisLength)* voxelSizeX;
            Objects.MedMinorAxisLenghtColoLCTN = median(ColoLCTNObjects.MinorAxisLength)* voxelSizeX;
    
            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoLCTNErodedMask = imerode(ColoLCTNMask, Conn6Strel);
            ColoLCTNPerimMask = (ColoLCTNMask - ColoLCTNErodedMask) > 0;

            % Additional feature images and vectors
            ColoLCTNBodyLabelIm = bwlabeln(ColoLCTNErodedMask, 6);

            % Erosion derived
            Objects.ColoLCTNPerimPixels = sum(ColoLCTNPerimMask(:));
            Objects.ColoLCTNBodyPixels = sum(ColoLCTNErodedMask(:)); 
            Objects.ColoLCTNBodyCount = max(ColoLCTNBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoLCTNShapeBySurface = Objects.ColoLCTNBodyPixels / Objects.ColoLCTNPerimPixels; % Roundness feature
            Objects.ColoLCTNBodycountByColoLCTNcount = Objects.ColoLCTNBodyCount / Objects.CountColoLCTN; % Invagination feature 
            % Previews
            ColoLCTNMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch4, [], 3)), imadjust(max(ch1, [], 3)));
            ColoLCTNMaskPreview_perimeter = imoverlay(ColoLCTNMaskPreview, bwperim(max(ColoLCTNMask, [], 3)), [1 1 1]);
            %imtool(ColoLCTNMaskPreview)
            SavePathColoLCTNperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoLCTN_peri.png'];
            imwrite(ColoLCTNMaskPreview_perimeter, SavePathColoLCTNperiPreview)
            else
            Objects.ColoLCTNArea = 0;
            Objects.ColoLCTNAreaNorm = 0;
            Objects.MinColoLCTNVol = 0;
            Objects.MaxColoLCTNVol = 0;
            Objects.MeanColoLCTNVol = 0;
            Objects.StdColoLCTNVol = 0;
            Objects.MedColoLCTNVol = 0;
            Objects.MadColoLCTNVol = 0;
            Objects.TotColoLCTNVolume = 0;
            Objects.TotColoLCTNVolumeNorm = 0;
            Objects.CountColoLCTN = 0;
            Objects.CountColoLCTNNorm = 0;                
			Objects.MedMajorAxisLenghtColoLCTN = 0;
            Objects.MedMinorAxisLenghtColoLCTN = 0;
    
            % Erosion derived
            Objects.ColoLCTNPerimPixels = 0;
            Objects.ColoLCTNBodyPixels = 0; 
            Objects.ColoLCTNBodyCount = 0; % Needed for invagination feature
            Objects.ColoLCTNShapeBySurface = 0; % Roundness feature
            Objects.ColoLCTNBodycountByColoLCTNcount = 0; % Invagination feature 
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
            ColoLATSMaskPreview = cat(3, imadjust(max(ch3, [], 3)), imadjust(max(ch4, [], 3)), imadjust(max(ch1, [], 3)));
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
            ColoLATNMaskPreview = cat(3, imadjust(max(ch3, [], 3)), imadjust(max(ch4, [], 3)), imadjust(max(ch1, [], 3)));
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

