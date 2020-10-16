function [Objects] = ImageAnalysisLC3_differentiation(ch1,ch2,InfoTableThis, PreviewPath)

         %% Area of cells
        Areaofcells = ch1>200; %it(Areaofcells)
        Areaofcells = imdilate(Areaofcells, strel('disk',3));
        %% Segmentation of all potential structures
        StructuresDoG = imfilter(ch2, fspecial('gaussian',50,1), 'symmetric') - imfilter(ch2, fspecial('gaussian',50,2), 'symmetric');
        StructuresMask = StructuresDoG>10;
        %% Segmentation of autophagic vacuoles -> similar approach to mitophagy since phagophores can't be recognized in the nuerons
        RatioIm = double(ch1) ./ double(ch2);
        StructuresCC = bwconncomp(StructuresMask, 26); % Connected components
        Objects = regionprops('table', StructuresCC, RatioIm, {'Area', 'PixelIdxList', 'PixelValues', 'MeanIntensity'});
        pHluorinThreshold = 1.1;
        Filtered_1 = Objects(Objects.Area < 200,:);
        ObjectsAutophagy = Filtered_1(Filtered_1.MeanIntensity < pHluorinThreshold, :);
        AutoPhagyMask = f_Create_Mask_from_ObjectList_Pixel_IDX(ObjectsAutophagy, 'PixelIdxList', ch1);
        AutoPhagyObjects = regionprops('table', bwconncomp(AutoPhagyMask), {'Area'}); 
        
        AreaofCells = regionprops('table', bwconncomp(Areaofcells), {'Area'});
        AreaofCells = sum(AreaofCells{:,1});
        % Add metadata
            BarCodeThis = unique(InfoTableThis.Barcode);
            AreaNameThis = InfoTableThis.AreaName;
            ColumnThis = InfoTableThis.Column;
            RowThis = InfoTableThis.ROW;
            FieldThis = InfoTableThis.field;
            TimepointThis= InfoTableThis.TimePoint;

            CellArrayDummy = cell(height(AutoPhagyObjects), 1);
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
            AreaofCellsColumn =  CellArrayDummy;
            AreaofCellsColumn = cellfun(@(x) AreaofCells, AreaofCellsColumn, 'UniformOutput', false);
            Objects = [AutoPhagyObjects, BarCodeColumn, AreaNameColumn, ColumnColumn, RowColumn, FieldColumn, TimePointColumn,AreaofCellsColumn];
            Objects.Properties.VariableNames(2:end) = {'Barcode','AreaName','Column','Row', 'Field','TimePoint','Areaofcells'};
            
            
            %% Previews
            
            imSize = [size(ch2,1),size(ch2,2)];
            [BarMask, ~] = f_barMask(20, 0.2152, imSize, imSize(1)-50, 50, 5); %it(BarMask)
            chEmpty = max(zeros(size(ch1),'uint16'),[],3);
            RGB = cat(3, imadjust(ch2), imadjust(ch1),  chEmpty); %it(RGB)

            ch12D = max(ch1, [], 3);
            ch22D = max(ch2, [], 3);
            
            PreviewDsRed = imoverlay(imadjust(ch22D), bwperim(max(AutoPhagyMask,[],3)), [1 0 0]);
            PreviewDsRed = imoverlay(PreviewDsRed, BarMask, [1 1 1]); %it(PreviewDsRed);
            
            PreviewPhluorin = imoverlay(imadjust(ch12D), bwperim(max(AutoPhagyMask,[],3)), [1 0 0]);
            PreviewPhluorin = imoverlay(PreviewPhluorin, BarMask, [1 1 1]); %it(PreviewPhluorin);
            
            PathRGB = [PreviewPath filesep BarCodeThis{:} '_t' num2str(TimepointThis) '_r' num2str(RowThis) '_c' num2str(ColumnThis) '_f' num2str(FieldThis) '_' AreaNameThis{:} '_RGB'];
            PathDsRed = [PreviewPath filesep BarCodeThis{:} '_t' num2str(TimepointThis) '_r' num2str(RowThis) '_c' num2str(ColumnThis) '_f' num2str(FieldThis) '_' AreaNameThis{:} '_ch2'];
            PathPHluorin = [PreviewPath filesep BarCodeThis{:} '_t' num2str(TimepointThis) '_r' num2str(RowThis) '_c' num2str(ColumnThis) '_f' num2str(FieldThis) '_' AreaNameThis{:} '_ch1'];
            

            imwrite(RGB, [PathRGB, '.png'])
            imwrite(PreviewDsRed, [PathDsRed, '.png'])
            imwrite(PreviewPhluorin, [PathPHluorin, '.png'])         
            
end

