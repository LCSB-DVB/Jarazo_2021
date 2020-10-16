function [] = f_CV8000_CreateSlideViewFiles(MetaData, SlideViewPlane, ThumbnailPath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Loop over slides (also corresponding to wells and timelines in this case)
    for TimeLine = 1:size(MetaData.TimeLines{:}, 1)

        %% Label organoids
        Xnorm = im2uint8(mat2gray(MetaData.TimeLines{1,1}.Fields{TimeLine,1}.X));
        Ynorm = im2uint8(mat2gray(MetaData.TimeLines{1,1}.Fields{TimeLine,1}.Y));

        OrganoidGroupIm = zeros(max(Ynorm)+1, max(Xnorm)+1, 'logical');
        for field = 1:length(Xnorm)
            OrganoidGroupIm((max(Ynorm)-Ynorm(field))+1, Xnorm(field)+1) = 1; % perspective from the microscope field display tool
        end

        %imtool(OrganoidGroupIm)

        LabelMatrix = bwlabeln(imclose(OrganoidGroupIm, strel('disk', 15)));
        %imtool(LabelMatrix,[])

        for field = 1:length(Xnorm)
            %MetaData.TimeLines{1,1}.Fields{TimeLine,1}(field, 'ROI_ID') = {LabelMatrix(Ynorm(field)+1, Xnorm(field)+1)};
            MetaData.TimeLines{1,1}.Fields{TimeLine,1}(field, 'ROI_ID') = {LabelMatrix((max(Ynorm)-Ynorm(field))+1, Xnorm(field)+1)};
        end

        %% Load organoid rescans
        InfoTable = MetaData.InfoTable{:};
        ImSize = size(imread(InfoTable{1,'file'}{:}));
        for o = 1:max(table2array(MetaData.TimeLines{1,1}.Fields{TimeLine,1}(:, 'ROI_ID'))) % organoids in this well (here timeline = well)
            FieldBoolThisOrganoid = table2array(MetaData.TimeLines{1,1}.Fields{TimeLine,1}(:, 'ROI_ID')) == o;
            FieldsThisOrganoid = sort(table2array(MetaData.TimeLines{1,1}.Fields{TimeLine,1}(FieldBoolThisOrganoid, 'Field')));
            CoordinatesThisOrganoid = MetaData.TimeLines{1,1}.Fields{TimeLine,1}(FieldBoolThisOrganoid, :);
            WellRowNum = MetaData.TimeLines{1,1}.Wells(TimeLine,:).WellRow;
            WellColNum = MetaData.TimeLines{1,1}.Wells(TimeLine,:).WellCol;
            ThisWellString = sprintf('%s%0.2d', char('A' + WellRowNum - 1), WellColNum);
            ImagePathsThisOrganoid = InfoTable(ismember(str2double(InfoTable.Field), CoordinatesThisOrganoid.Field) & strcmp(InfoTable.Well,ThisWellString),:);
            FieldCoordinates = MetaData.TimeLines{1,1}.Fields{TimeLine,1};
            XVecThisOrganoid = sort(unique(table2array(FieldCoordinates(ismember(FieldCoordinates.Field, FieldsThisOrganoid), 'X'))));
            YVecThisOrganoid = sort(unique(table2array(FieldCoordinates(ismember(FieldCoordinates.Field, FieldsThisOrganoid), 'Y'))));
            Im = cell(4,1);
            for ch = 1:4
                Im{ch} = zeros(length(YVecThisOrganoid) * ImSize(1), length(XVecThisOrganoid) * ImSize(2));
                for field = 1:length(FieldsThisOrganoid)
                    FieldThis = FieldsThisOrganoid(length(FieldsThisOrganoid)+1 - field); % start from last acquired field
                    %ThisTileInfo = ImagePathsThisOrganoid(str2double(ImagePathsThisOrganoid.Channel) == ch & str2double(ImagePathsThisOrganoid.Field) == FieldThis, :);
                    ThisTileInfo = ImagePathsThisOrganoid(str2double(ImagePathsThisOrganoid.Channel) == ch & ...
                                                          str2double(ImagePathsThisOrganoid.Field) == FieldThis & ...
                                                          strcmp(ImagePathsThisOrganoid.Plane, num2str(SlideViewPlane)), :);
                    ThisTilePath = ThisTileInfo{1,'file'}{:};
                    TileThis = imread(ThisTilePath);
                    %StartRowThisTile =    ((length(YVecThisOrganoid) - (find(YVecThisOrganoid == table2array(FieldCoordinates(FieldThis, 'Y'))) - 1)) * ImSize(1)) + 1;
                    StartRowThisTile =    ((length(YVecThisOrganoid) - (find(YVecThisOrganoid == table2array(FieldCoordinates(FieldThis, 'Y'))))) * ImSize(1)) + 1;
                    StartColumnThisTile = ((find(XVecThisOrganoid == table2array(FieldCoordinates(FieldThis, 'X'))) - 1) * ImSize(2)) + 1;
                    Im{ch}(StartRowThisTile:StartRowThisTile+ImSize(1)-1, StartColumnThisTile:StartColumnThisTile+ImSize(2)-1) = TileThis;
                    %% Save compact Organoid channel image including top left pixel coordinate and slide size to file
                    OrganoidThumbNameThis = sprintf('%s/%s_organoidID%02d_channel%d.mat', ThumbnailPath, ThisWellString, o, ch);
                    ShrinkFactor = 20;
                    OrganoidImThis = imresize(Im{ch}, 1/ShrinkFactor);% imtool(OrganoidImThis, [])
                    save(OrganoidThumbNameThis, 'OrganoidImThis', 'OrganoidGroupIm', 'LabelMatrix')

                end
    %             if TimeLine == 2 & o == 1
    %                 imtool(Im{ch},[])
    %             end
            end
        end
    end % TimeLine
    
    disp('Successfully created files for slideview.')

end

