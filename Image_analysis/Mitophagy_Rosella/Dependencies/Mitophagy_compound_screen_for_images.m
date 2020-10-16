function [Objects] = Mitophagy_compound_screen(files, index, InfoTable, PreviewPath)
%Mitophagy analysis Javier Jonathan
%   Detailed explanation goes here


            cube = readflexcube(files{index}, 'PlaneCount', 5); % Read 4-D image cube            
            ch1 = cube.data(:, :, :, 1); %pHLourin vol(ch1, 0, 1000)
            ch2 = cube.data(:, :, :, 2); %DsRed vol(ch2, 0, 1000)
            %ch3 = cube.data(:, :, :, 3); % vol(ch3, 0, 100)
 %% Raw Images
 % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(20, 0.2152, imSize, imSize(1)-50, 550, 10);
    %it(BarMask)
    
    %Im_pHluorinImRaw = imadjust(max(ch1, [], 3), [0 0.02], [0 1]); 
    %imtool(Im_pHluorinImRaw); %imwrite(Im_pHluorinImRaw,'Im_pHluorinImRaw.png')
    %[ZoomOverlay Zoom] = f_ZoomBoxRGB(Im_pHluorinImRaw, 180, 240, 50, 50, 3, 1);
    % ZoomOverlay = imoverlay2(ZoomOverlay, BarMask, [1 1 1]);
    %it(ZoomOverlay) %imwrite(ZoomOverlay,'Im_pHluorinImRaw_ZO.png')
    %it(Zoom) %imwrite(Zoom,'Im_pHluorinImRaw_Z.png')
    
    %Im_DsRedImRaw = imadjust(max(ch2, [], 3), [0 0.009], [0 1]); 
    %imtool(Im_DsRedImRaw); %imwrite(Im_DsRedImRaw,'Im_DsRedImRaw.png')
    %[ZoomOverlay Zoom] = f_ZoomBoxRGB(Im_DsRedImRaw, 180, 240, 50, 50, 3, 1);
    % ZoomOverlay = imoverlay2(ZoomOverlay, BarMask, [1 1 1]);
    %it(ZoomOverlay) %imwrite(ZoomOverlay,'Im_DsRedImRaw_ZO.png')
    %it(Zoom) %imwrite(Zoom,'Im_DsRedImRaw_Z.png')          
            
    %Segment mitochondria
            MitoDoGa = imfilter(ch2, fspecial('gaussian',50,1), 'symmetric') - imfilter(ch2, fspecial('gaussian',50,2), 'symmetric');
            %Im_MitoDoG = imadjust(max(MitoDoGa, [], 3), [0 0.0005], [0 1]); 
            %imtool(Im_MitoDoG); %imwrite(Im_MitoDoG,'Im_MitoDoG.png')
            %[ZoomOverlay Zoom] = f_ZoomBoxRGB(Im_MitoDoG, 180, 240, 50, 50, 3, 1);
            %it(ZoomOverlay) %imwrite(ZoomOverlay,'Im_MitoDoG_ZO.png')
            %it(Zoom) %imwrite(Zoom,'Im_MitoDoG_Z.png')
            
            MitoDoGb = imfilter(ch2, fspecial('gaussian',50,1), 'symmetric') - imfilter(ch2, fspecial('gaussian',50,5), 'symmetric');
            %vol(MitoDoGb,0, 50)
            %vol(MitoDoGa, 0, 10)
            %vol(MitoDoGb)
            
            MitoMask = MitoDoGa > 10; % vol(MitoMask)
            %Im_MitoMask = max(MitoMask, [], 3); 
            %imtool(Im_MitoMask); %imwrite(Im_MitoMask,'Im_MitoMask.png')
            %[ZoomOverlay Zoom] = f_ZoomBoxRGB(Im_MitoMask, 180, 240, 50, 50, 3, 1);
            %it(ZoomOverlay) %imwrite(ZoomOverlay,'Im_MitoMask_ZO.png')
            %it(Zoom) %imwrite(Zoom,'Im_MitoMask_Z.png')

            RatioIm = double(ch1) ./ double(ch2); %vol(RatioIm)
            %Im_RatioIm = max(RatioIm, [], 3); 
            %Manually adjust Window of values Min:-0.1081 Max:4.7391
            % Save the picture manually
            % load the picture -> Im_RatioIm2= imread('Im_RatioIm2.png');
            %imtool(Im_RatioIm2); %imwrite(Im_RatioIm2,'Im_RatioIm.png')
            %[ZoomOverlay Zoom] = f_ZoomBoxRGB(Im_RatioIm2, 180, 240, 50, 50, 3, 1);
            %imtool(ZoomOverlay) %imwrite(ZoomOverlay,'Im_RatioIm_ZO.png')
            %it(Zoom) %imwrite(Zoom,'Im_RatioIm_Z.png')
            
            MitoCC = bwconncomp(MitoMask, 26); % Connected components
                 
            MitoObjects = regionprops('table', MitoCC, RatioIm, {'Area', 'PixelIdxList', 'PixelValues', 'MeanIntensity'});
            
            % figure; histogram(MitoObjects.MeanIntensity, 100)
            % pHluorinThreshold = 0.6;
            pHluorinThreshold = 1;
            MitoObjectsMitophagy = MitoObjects(MitoObjects.MeanIntensity < pHluorinThreshold, :);
            MitoObjectsBaseLevel = MitoObjects(MitoObjects.MeanIntensity >= pHluorinThreshold, :);
            MitoPhagyMask = f_Create_Mask_from_ObjectList_Pixel_IDX(MitoObjectsMitophagy, 'PixelIdxList', ch1);
            MitoBaseMask = f_Create_Mask_from_ObjectList_Pixel_IDX(MitoObjectsBaseLevel, 'PixelIdxList', ch1) * 2;
            
            
            % vol(MitoPhagyMask + MitoBaseMask, 0,2,'jet')
            % vol(ch3,0,100)
            % vol(ch2,0,100)
            %vol(MitoDoGb, 0,100)
            
            %MitophagyReconstructedMask = imreconstruct(logical(MitoPhagyMask), logical(MitoDoGb > 50)); %vol(MitophagyReconstructedMask)
            MitophagyReconstructedMask = logical(MitoPhagyMask);
            
            MitoBaseMask = bwareaopen(MitoBaseMask, 2); %vol(MitoBaseMask)
            %Im_Mitochondria = max(MitoBaseMask, [], 3); 
            %imtool(Im_Mitochondria); %imwrite(Im_Mitochondria,'Im_Mitochondria.png')
            %[ZoomOverlay Zoom] = f_ZoomBoxRGB(Im_Mitochondria, 180, 240, 50, 50, 3, 1);
            %it(ZoomOverlay) %imwrite(ZoomOverlay,'Im_Mitochondria_ZO.png')
            %it(Zoom) %imwrite(Zoom,'Im_Mitochondria_Z.png')
            
            MitophagyReconstructedMask = bwareaopen(MitophagyReconstructedMask, 2);
            %Im_Mitophagy_Event = max(MitophagyReconstructedMask, [], 3); 
            %imtool(Im_Mitophagy_Event); %imwrite(Im_Mitophagy_Event,'Im_Mitophagy Event.png')
            %[ZoomOverlay Zoom] = f_ZoomBoxRGB(Im_Mitophagy_Event, 180, 240, 50, 50, 3, 1);
            %it(ZoomOverlay) %imwrite(ZoomOverlay,'Im_Mitophagy_Event_ZO.png')
            %it(Zoom) %imwrite(Zoom,'Im_Mitophagy_Event_Z.png')
            
            %% Collect Features
            MitoObjects = regionprops('table', bwconncomp(MitoBaseMask), {'Area'}); % pixel count in 3D is volume
            MitoClass = cell(height(MitoObjects),1); MitoClass = cellfun(@(x) 1, MitoClass);
            MitoObjects = [MitoObjects, array2table(MitoClass)];
            MitoPhagyObjects = regionprops('table', bwconncomp(MitophagyReconstructedMask), {'Area'}); % pixel count in 3D is volume
            MitoClass = cell(height(MitoPhagyObjects),1); MitoClass = cellfun(@(x) 2, MitoClass);
            MitoPhagyObjects = [MitoPhagyObjects, array2table(MitoClass)];
            Objects = [MitoObjects; MitoPhagyObjects];
            
            % Add metadata
            BarCodeThis = unique(InfoTable.Barcode);
            AreaNameThis = InfoTable.AreaName(index);
            ColumnThis = InfoTable.Column(index);
            RowThis = InfoTable.Row(index);
            FieldThis = InfoTable.field(index);
            TimepointThis= InfoTable.TimePoint(index);

            CellArrayDummy = cell(height(Objects), 1);
            BarCodeColumn = CellArrayDummy;
            BarCodeColumn = cellfun(@(x) BarCodeThis, BarCodeColumn, 'UniformOutput', false);
            AreaNameColumn = CellArrayDummy;
            AreaNameColumn = cellfun(@(x) AreaNameThis, AreaNameColumn, 'UniformOutput', false);        
            ColumnColumn = CellArrayDummy;
            ColumnColumn = cellfun(@(x) ColumnThis, ColumnColumn, 'UniformOutput', false);
            RowColumn = CellArrayDummy;
            RowColumn = cellfun(@(x) RowThis, RowColumn, 'UniformOutput', false);
            FieldColumn = CellArrayDummy;
            FieldColumn = cellfun(@(x) FieldThis, FieldColumn, 'UniformOutput', false);
            TimePointColumn = CellArrayDummy;
            TimePointColumn = cellfun(@(x) TimepointThis, TimePointColumn, 'UniformOutput', false);
            
            Objects = [Objects, BarCodeColumn, AreaNameColumn, ColumnColumn, RowColumn, FieldColumn, TimePointColumn];
            Objects.Properties.VariableNames(3:end) = {'Barcode','AreaName','Column','Row', 'Field','TimePoint'};
            
            %% Previews
            
            imSize = [size(ch2,1),size(ch2,2)];
            [BarMask, BarCenter] = f_barMask(20, 0.2152, imSize, imSize(1)-50, 50, 5); %it(BarMask)
            chEmpty = max(zeros(size(ch1),'uint16'),[],3);
            RGB = cat(3, imadjust(ch2(:,:,5)), imadjust(ch1(:,:,5)),  chEmpty); %it(RGB)

            ch12D = max(ch1, [], 3);
            ch22D = max(ch2, [], 3);
            
            MitoMask2D = max(MitoBaseMask + MitophagyReconstructedMask, [], 3);

            PreviewDsRed = imoverlay(imadjust(ch22D), bwperim(max(MitoBaseMask,[],3)), [0 1 0]);
            PreviewDsRed = imoverlay(PreviewDsRed, bwperim(max(MitophagyReconstructedMask,[],3)), [1 0 0]);
            PreviewDsRed = imoverlay(PreviewDsRed, BarMask, [1 1 1]); %it(PreviewDsRed);
            
            PreviewPhluorin = imoverlay(imadjust(ch12D), bwperim(max(MitoBaseMask,[],3)), [0 1 0]);
            PreviewPhluorin = imoverlay(PreviewPhluorin, bwperim(max(MitophagyReconstructedMask,[],3)), [1 0 0]);
            PreviewPhluorin = imoverlay(PreviewPhluorin, BarMask, [1 1 1]); %it(PreviewPhluorin);
            %%
            PathRGB = [PreviewPath filesep BarCodeThis{:} '_t' num2str(TimepointThis) '_r' num2str(RowThis) '_c' num2str(ColumnThis) '_f' num2str(FieldThis) '_' AreaNameThis{:} '_RGB'];
            PathDsRed = [PreviewPath filesep BarCodeThis{:} '_t' num2str(TimepointThis) '_r' num2str(RowThis) '_c' num2str(ColumnThis) '_f' num2str(FieldThis) '_' AreaNameThis{:} '_ch2'];
            PathPHluorin = [PreviewPath filesep BarCodeThis{:} '_t' num2str(TimepointThis) '_r' num2str(RowThis) '_c' num2str(ColumnThis) '_f' num2str(FieldThis) '_' AreaNameThis{:} '_ch1'];
            

            imwrite(RGB, [PathRGB, '.png'])
            imwrite(PreviewDsRed, [PathDsRed, '.png'])
            imwrite(PreviewPhluorin, [PathPHluorin, '.png'])
            
%             %% More previews for the workflow figure
%             PreviewMitoDoG = imadjust(max(MitoDoGa, [], 3), [0 0.005], [0 1]);
%             PreviewMitoDoG = imoverlay(PreviewMitoDoG, BarMask, [1 1 1]);
%             %imtool(PreviewMitoDoG)
%             
%             PreviewMitoMask = imoverlay(imadjust(MitoMask2D), BarMask, [1 1 1]);
%             %imtool(PreviewMitoMask)
%             
%             PreviewMitophagySeedMask = imoverlay(imadjust(max(MitoPhagyMask, [], 3)), BarMask, [1 1 1]);
%             %imtool(PreviewMitophagySeedMask)
%             
%             PreviewMitophagyLimitingMask = imoverlay(max(MitophagyReconstructedMask, [], 3), BarMask, [1 1 1]);
%             %imtool(PreviewMitophagyLimitingMask)
%             
%             PathMitoDoG = [PreviewPath filesep BarCodeThis{:} '_' num2str(RowThis) '_' num2str(ColumnThis) '_' num2str(FieldThis) '_' AreaNameThis{:} '_MitoDoG'];
%             PathMitoMask = [PreviewPath filesep BarCodeThis{:} '_' num2str(RowThis) '_' num2str(ColumnThis) '_' num2str(FieldThis) '_' AreaNameThis{:} '_MitoMask'];
%             PathMitophagySeedMask = [PreviewPath filesep BarCodeThis{:} '_' num2str(RowThis) '_' num2str(ColumnThis) '_' num2str(FieldThis) '_' AreaNameThis{:} '_MitophagySeedMask'];
%             PathMitophagyLimitingMask = [PreviewPath filesep BarCodeThis{:} '_' num2str(RowThis) '_' num2str(ColumnThis) '_' num2str(FieldThis) '_' AreaNameThis{:} '_MitophagyLimitingMask'];
%             
%             imwrite(PreviewMitoDoG, [PathMitoDoG, '.png'])
%             imwrite(PreviewMitoMask, [PathMitoMask, '.png'])
%             imwrite(PreviewMitophagySeedMask, [PathMitophagySeedMask, '.png'])
%             imwrite( PreviewMitophagyLimitingMask, [PathMitophagyLimitingMask, '.png'])
            
            
        
end

