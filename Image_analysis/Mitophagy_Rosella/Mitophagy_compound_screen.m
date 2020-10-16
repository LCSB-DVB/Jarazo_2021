function [Objects] = Mitophagy_compound_screen(files, index, InfoTable, PreviewPath)

%   Detailed explanation goes here


            cube = readflexcube(files{index}, 'PlaneCount', 5); % Read 4-D image cube            
            ch1 = cube.data(:, :, :, 1); %pHLourin vol(ch1, 0, 1000)
            ch2 = cube.data(:, :, :, 2); %DsRed vol(ch2, 0, 1000)
            
            
            %Segment mitochondria
            MitoDoGa = imfilter(ch2, fspecial('gaussian',50,1), 'symmetric') - imfilter(ch2, fspecial('gaussian',50,2), 'symmetric');
            %vol(MitoDoGa, 0, 10)
            %vol(MitoDoGb)
            
            MitoMask = MitoDoGa > 10; % vol(MitoMask)


            RatioIm = double(ch1) ./ double(ch2); %vol(RatioIm)
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
            
            MitophagyReconstructedMask = logical(MitoPhagyMask);
            MitoBaseMask = bwareaopen(MitoBaseMask, 2); %vol(MitoBaseMask)
            MitophagyReconstructedMask = bwareaopen(MitophagyReconstructedMask, 2);
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
            [BarMask, ~] = f_barMask(20, 0.2152, imSize, imSize(1)-50, 50, 5); %it(BarMask)
            chEmpty = max(zeros(size(ch1),'uint16'),[],3);
            RGB = cat(3, imadjust(ch2(:,:,5)), imadjust(ch1(:,:,5)),  chEmpty); %it(RGB)

            ch12D = max(ch1, [], 3);
            ch22D = max(ch2, [], 3);
            
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

