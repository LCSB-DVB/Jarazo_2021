function [Objects] = Analysis_TFEB_neurons(InfoTableThis,  PreviewPath, ch1, ch2, ch3)

   %% image analysis
%it(ch1) Nuclei
%it(ch2) TFEB
%it(ch3) TH
%% Segment nuclei
             
	NucleiMask = ch1 > 200;
    NucleiMask = bwareaopen(NucleiMask, 10); %it(NucleiMask); it(ch1);
        if sum(NucleiMask(:))== 0
            Objects = {}; 
            return
        end
    NucleiDead = ch1 > 800; 
    NucleiDead = imdilate (NucleiDead, strel('sphere',2));%it(NucleiDead)
    NucleiLive = NucleiMask & ~NucleiDead; %it(NucleiLive)
    NucleiLive = bwareaopen(NucleiLive, 10);
    
%% Segment TFEB
    TFEBBlurred = imfilter(ch2, fspecial('gaussian',5, 4)); % it(TFEBBlurred)    
    TFEBMask = TFEBBlurred > 200;
    TFEBMask = bwareaopen(TFEBMask, 20);% it(TFEBMask)
	
%% Segment TH soma
    TH_SomaDoG = imfilter(ch3, fspecial('gaussian', 31, 8) - fspecial('gaussian', 31, 10), 'symmetric'); %it(TH_SomaDoG)
    
    THSomaMask = TH_SomaDoG > 35;  % it(THSomaMask)
    THSomaMask = bwareaopen(THSomaMask, 50);
    THSomaMask = imdilate (THSomaMask, strel('sphere',2));


%% Colocalization Nuclei Mask
    ColoNucMask = NucleiMask & TFEBMask; %it(ColoNucMask)   
    
%% Colocalization TH-TFEB
	ColoTHMask = THSomaMask & ColoNucMask; %it(ColoTHMask)    

%% Morphometrics
    TFEBLM = bwlabeln(TFEBMask);
    TFEBObjects = regionprops('table', TFEBLM, {'Area'});
    NucleiLM = bwlabeln(NucleiMask);
    NucleiObjects = regionprops('table', NucleiLM, {'Area'});
    ColoNucLM = bwlabeln(ColoNucMask);
    ColoNucObjects = regionprops('table', ColoNucLM, {'Area'});
	ColoTHLM = bwlabeln(ColoTHMask);
    ColoTHObjects = regionprops('table', ColoTHLM, {'Area'});                        

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
    Objects.NucAliveNorm = sum(NucleiLive(:))/sum(NucleiMask(:)); 
    Objects.NucDeadNorm = sum(NucleiDead(:))/sum(NucleiMask(:)); 
	%% TFEB derived
    Objects.TFEBArea = sum(TFEBMask(:));
    Objects.TFEBAreaNorm = sum(TFEBMask(:))/sum(NucleiMask(:));
    Objects.CountTFEBObj = size(TFEBObjects,1);
	Objects.CountTFEBObjNorm = size(TFEBObjects,1)/size(NucleiObjects,1);
    TFEBInten = ch2 .* uint16(TFEBMask);
    Objects.TFEBIntenNorm = sum(TFEBInten(:))/sum(NucleiMask(:));    

%% 2D previews
    
    % Scalebar
    imSize = size(ch1);
    voxelSizeX = 0.646;
	[BarMask, ~] = f_barMask(20, voxelSizeX, imSize, imSize(1)-50, 75, 10);
    %it(BarMask)
    RGB = cat(3, imadjust(ch2), imadjust(ch3), imadjust(ch1));
    RGB = imoverlay(RGB, BarMask, [1 1 1]);
    % imtool(RGB)
  
    TFEBMaskPreview = imoverlay(imadjust(ch2), bwperim(TFEBMask), [1 0 0]);
    TFEBMaskPreview = imoverlay(TFEBMaskPreview, BarMask, [1 1 1]);
    %imtool(TFEBMaskPreview )
    TFEBRawPreview = imadjust(ch2);
    TFEBRawPreview = imoverlay(TFEBRawPreview, BarMask, [1 1 1]);
    %imtool(TFEBRawPreview)
    
    
    SavePathTFEBMaskPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_TFEB_mask.png'];
    SavePathTFEBRawPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_TFEB_raw.png'];
    SavePathRGBPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RGB.png'];
    
   
    imwrite(TFEBMaskPreview, SavePathTFEBMaskPreview)
    imwrite(TFEBRawPreview,  SavePathTFEBRawPreview)
    imwrite(RGB, SavePathRGBPreview)
 
 
 %% Colo derived
    if sum(ColoNucMask(:))> 0 && sum(ColoTHMask(:))> 0
				% ColoNuc derived
				Objects.ColoNucMaskArea = sum(ColoNucMask(:));
				Objects.ColoNucMaskNorm = sum(ColoNucMask(:))/sum(NucleiMask(:));
				Objects.CountColoNucObj = size(ColoNucObjects,1);
				Objects.CountColoNucObjNorm = size(ColoNucObjects,1)/size(NucleiObjects,1);
				ColoNucInten = ch2 .* uint16(ColoNucMask);
				Objects.ColoNucIntenNorm = sum(ColoNucInten(:))/sum(NucleiMask(:));
				% ColoTH derived
				Objects.ColoTHMaskArea = sum(ColoTHMask(:));
				Objects.ColoTHMaskNorm = sum(ColoTHMask(:))/sum(NucleiMask(:));
				Objects.CountColoTHObj = size(ColoTHObjects,1);
				Objects.CountColoTHObjNorm = size(ColoTHObjects,1)/size(NucleiObjects,1);
				ColoTHInten = ch2 .* uint16(ColoTHMask);
				Objects.ColoTHIntenNorm = sum(ColoTHInten(:))/sum(NucleiMask(:));
							
				
				% Previews
				ColoNucMaskPeriPreview = imoverlay(imadjust(ch2), bwperim(ColoNucMask), [1 0 0]);
				ColoNucMaskPeriPreview = imoverlay(ColoNucMaskPeriPreview, BarMask, [1 1 1]);
				%imtool(ColoNucMaskPeriPreview )
				ColoNucMaskPreview = imoverlay(imadjust(ch2), ColoNucMask, [0 1 0]);
				ColoNucMaskPreview = imoverlay(ColoNucMaskPreview, BarMask, [1 1 1]);
				%imtool(ColoNucMaskPreview )
				
				ColoTHMaskPeriPreview = imoverlay(imadjust(ch2), bwperim(ColoTHMask), [1 0 0]);
				ColoTHMaskPeriPreview = imoverlay(ColoTHMaskPeriPreview, BarMask, [1 1 1]);
				%imtool(ColoTHMaskPeriPreview )
				ColoTHMaskPreview = imoverlay(imadjust(ch2), ColoTHMask, [0 1 0]);
				ColoTHMaskPreview = imoverlay(ColoTHMaskPreview, BarMask, [1 1 1]);
				%imtool(ColoTHMaskPreview )
				
				SavePathColoNucPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoNuc.png'];
				SavePathColoNucperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoNuc_peri.png'];
				
				imwrite(ColoNucMaskPreview, SavePathColoNucPreview)
				imwrite(ColoNucMaskPeriPreview, SavePathColoNucperiPreview)
				
				SavePathColoTHPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoTH.png'];
				SavePathColoTHperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoTH_peri.png'];
				
				imwrite(ColoTHMaskPreview, SavePathColoTHPreview)
				imwrite(ColoTHMaskPeriPreview, SavePathColoTHperiPreview)
				
    elseif sum(ColoNucMask(:))> 0 && sum(ColoTHMask(:))==0
				% ColoNuc derived
				Objects.ColoNucMaskArea = sum(ColoNucMask(:));
				Objects.ColoNucMaskNorm = sum(ColoNucMask(:))/sum(NucleiMask(:));
				Objects.CountColoNucObj = size(ColoNucObjects,1);
				Objects.CountColoNucObjNorm = size(ColoNucObjects,1)/size(NucleiObjects,1);
				ColoNucInten = ch2 .* uint16(ColoNucMask);
				Objects.ColoNucIntenNorm = sum(ColoNucInten(:))/sum(NucleiMask(:));
				% ColoTH derived
				Objects.ColoTHMaskArea = 0;
				Objects.ColoTHMaskNorm = 0;
				Objects.CountColoTHObj = 0;
				Objects.CountColoTHObjNorm = 0;
				Objects.ColoTHIntenNorm = 0;
							
				
				% Previews
				ColoNucMaskPeriPreview = imoverlay(imadjust(ch2), bwperim(ColoNucMask), [1 0 0]);
				ColoNucMaskPeriPreview = imoverlay(ColoNucMaskPeriPreview, BarMask, [1 1 1]);
				%imtool(ColoNucMaskPeriPreview )
				ColoNucMaskPreview = imoverlay(imadjust(ch2), ColoNucMask, [0 1 0]);
				ColoNucMaskPreview = imoverlay(ColoNucMaskPreview, BarMask, [1 1 1]);
				%imtool(ColoNucMaskPreview )
										
				SavePathColoNucPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoNuc.png'];
				SavePathColoNucperiPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_ColoNuc_peri.png'];
				
				imwrite(ColoNucMaskPreview, SavePathColoNucPreview)
				imwrite(ColoNucMaskPeriPreview, SavePathColoNucperiPreview)
            
	else
				% ColoNuc derived
				Objects.ColoNucMaskArea = 0;
				Objects.ColoNucMaskNorm = 0;
				Objects.CountColoNucObj = 0;
				Objects.CountColoNucObjNorm = 0;
				Objects.ColoNucIntenNorm = 0;
				% ColoTH derived
				Objects.ColoTHMaskArea = 0;
				Objects.ColoTHMaskNorm = 0;
				Objects.CountColoTHObj = 0;
				Objects.CountColoTHObjNorm = 0;
				Objects.ColoTHIntenNorm = 0;
	
    end    
end

