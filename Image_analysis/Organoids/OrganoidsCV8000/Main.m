display(['directory: ', pwd])
%% User inputs
%SelectOrganoids = 'Unknown';
SelectOrganoids = 'Known';
PlaneCount = 25;
delete(gcp('nocreate'))
pool = parpool(28);
%addpath(genpath('/mnt/lscratch/users/pantony/GitLibraries'))
addpath(genpath('/work/projects/lcsb_hcs/Library/hcsforge'))
InPath = '/mnt/lscratch/users/jjarazo/HCSdata/JJ_20190729_PINK1_N_A1-A4_20190729_111528/JJ_20190729_PINK1_slide_1-4_exp4-5_07-29-19_11-22-26';
DataPath = InPath;
OutPath = '/mnt/lscratch/users/jjarazo/Results/PINK1_slide_1-4_exp4-5/Outputs';

%% Prepare folders
mkdir(OutPath)
ThumbnailPath = [OutPath, filesep, 'Thumbnails'];
mkdir(ThumbnailPath)
PreviewPath = [OutPath, filesep, 'Previews'];
mkdir(PreviewPath)


%% Load Metadata
mesFile = dir('./*mes');
mesFile = mesFile.name;
MetaData = f_CV8000_getChannelInfo(DataPath, mesFile);

if strcmp(SelectOrganoids, 'Unknown')
    tic
    disp('Creating slide preview to assist you in choosing organoids and sample naming.')
    SlideViewPlane = 12;
    f_CV8000_CreateSlideViewFiles(MetaData, SlideViewPlane, ThumbnailPath)

    %% Reconstruct slide A4 channel 1
    GrayRangeInput = [0 0.03]; % See imadjust range [0 1]
    SlidePreviewIms = f_CV8000_SlideViewOrganoids(ThumbnailPath, unique(MetaData.InfoTable{:}.Well), 1, GrayRangeInput);
    for s = 1:size(SlidePreviewIms, 1)
        %imtool(SlidePreviewIms{s})% Index corresponds to slide
        imwrite(SlidePreviewIms{s}, [ThumbnailPath,filesep,'SlidePreview_',num2str(s),'.png'])
    end
    toc
%     waitfor(msgbox('Please create a txt file in the working directory that is space delimited, where the first column is the Matlab organoid ID as shown in this slideview, the second column the Well in format A01 and the third column is the sample name. The name of this file needs to be assigned to the variable SelectOrganoids.'))
%     disp('Please trigger the script again after updating SelectOrganoids')
    return
end


%% Analyze selected organoids
InfoTable = MetaData.InfoTable{:};
OrganoidsToAnalyze = readtable('./OrganoidSelection.txt','Delimiter',' ','ReadVariableNames',false);
OrganoidsToAnalyze.Properties.VariableNames = {'Idx', 'Well', 'AreaName'};
OrganoidsToAnalyze.OrganoidID = rowfun(@(a,b) {sprintf('%s_%02d', a{:}, b)}, OrganoidsToAnalyze, 'InputVariables', {'Well', 'Idx'});

%% Identify correspondance between wells and timelines
TimeLines = MetaData.TimeLines{:};
TimeLineWells = {};
for t = 1:size(TimeLines, 1)
    TimeLineWells{t} = sprintf('%s%02d', char('A'+TimeLines{t, 'Wells'}.WellRow-1), TimeLines{t, 'Wells'}.WellCol);
end

ObjectsAll = {};
for organoid = 1:size(OrganoidsToAnalyze,1)
%for organoid = 1:2
    tic
    
    OrganoidThis = OrganoidsToAnalyze(organoid, :);
    SlideThis = OrganoidThis.Well{:};
    Sample = OrganoidThis.AreaName{:};
    ID = OrganoidThis.Idx;
    disp(['Well ', OrganoidsToAnalyze{organoid, 'Well'}{:}, '__Organoid #' num2str(organoid)])
    ThumbnailPathThis = [ThumbnailPath, filesep, SlideThis, '_organoidID', sprintf('%02d',ID), '_channel', num2str(1), '.mat']; % Defaults to channel 1    
    ThumbnailThis = load(ThumbnailPathThis);
    
    %% Find timeline corresponding to this slide
    TimeLine = find(strcmp(TimeLineWells, OrganoidsToAnalyze{organoid, 'Well'}));

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
        MetaData.TimeLines{1,1}.Fields{TimeLine,1}(field, 'ROI_ID') = {LabelMatrix((max(Ynorm)-Ynorm(field))+1, Xnorm(field)+1)};
    end

    %% Load organoid rescan image mosaic (optimized for speed)
    ImSize = size(imread(InfoTable{1,'file'}{:}));
    FieldBoolThisOrganoid = table2array(MetaData.TimeLines{1,1}.Fields{TimeLine,1}(:, 'ROI_ID')) == ID;
    FieldsThisOrganoid = sort(table2array(MetaData.TimeLines{1,1}.Fields{TimeLine,1}(FieldBoolThisOrganoid, 'Field')));
    CoordinatesThisOrganoid = MetaData.TimeLines{1,1}.Fields{TimeLine,1}(FieldBoolThisOrganoid, :);
    WellRowNum = MetaData.TimeLines{1,1}.Wells(TimeLine,:).WellRow;
    WellColNum = MetaData.TimeLines{1,1}.Wells(TimeLine,:).WellCol;
    ThisWellString = sprintf('%s%0.2d', char('A' + WellRowNum - 1), WellColNum);
    ImagePathsThisOrganoid = InfoTable(ismember(str2double(InfoTable.Field), CoordinatesThisOrganoid.Field) & strcmp(InfoTable.Well,ThisWellString),:);
    FieldCoordinates = MetaData.TimeLines{1,1}.Fields{TimeLine,1};
    XVecThisOrganoid = sort(unique(table2array(FieldCoordinates(ismember(FieldCoordinates.Field, FieldsThisOrganoid), 'X'))));
    YVecThisOrganoid = sort(unique(table2array(FieldCoordinates(ismember(FieldCoordinates.Field, FieldsThisOrganoid), 'Y'))));
    ch1Meta=table();
    ch2Meta=table();
    ch3Meta=table();
    ch4Meta=table();    
    progressCh1 = 0;
    progressCh2 = 0;
    progressCh3 = 0;
    progressCh4 = 0;
    
    for ch = 1:4
        for field = 1:length(FieldsThisOrganoid)
            for plane = 1:PlaneCount
                FieldThis = FieldsThisOrganoid(length(FieldsThisOrganoid)+1 - field); % start from last acquired field
                ThisTileInfo = ImagePathsThisOrganoid(str2double(ImagePathsThisOrganoid.Channel) == ch & ...
                                                      str2double(ImagePathsThisOrganoid.Field) == FieldThis & ...
                                                      strcmp(ImagePathsThisOrganoid.Plane, sprintf('%02d', plane)), :);
                ThisTilePath = ThisTileInfo{1,'file'}{:};
                StartRowThisTile =    ((length(YVecThisOrganoid) - (find(YVecThisOrganoid == table2array(FieldCoordinates(FieldThis, 'Y'))))) * ImSize(1)) + 1;
                StartColumnThisTile = ((find(XVecThisOrganoid == table2array(FieldCoordinates(FieldThis, 'X'))) - 1) * ImSize(2)) + 1;
                switch ch
                    case 1
                        progressCh1 = progressCh1 + 1;
                        ch1Meta(progressCh1, 'StartRowThisTile') = {StartRowThisTile};
                        ch1Meta(progressCh1, 'StartColumnThisTile') = {StartColumnThisTile};
                        ch1Meta(progressCh1, 'Plane') = {plane};
                        ch1Meta(progressCh1, 'Field') = {field};
                        ch1Meta(progressCh1, 'ThisTilePath') = {ThisTilePath};
                    case 2
                        progressCh2 = progressCh2 + 1;
                        ch2Meta(progressCh2, 'StartRowThisTile') = {StartRowThisTile};
                        ch2Meta(progressCh2, 'StartColumnThisTile') = {StartColumnThisTile};
                        ch2Meta(progressCh2, 'Plane') = {plane};
                        ch2Meta(progressCh2, 'Field') = {field};
                        ch2Meta(progressCh2, 'ThisTilePath') = {ThisTilePath};
                    case 3
                        progressCh3 = progressCh3 + 1;
                        ch3Meta(progressCh3, 'StartRowThisTile') = {StartRowThisTile};
                        ch3Meta(progressCh3, 'StartColumnThisTile') = {StartColumnThisTile};
                        ch3Meta(progressCh3, 'Plane') = {plane};
                        ch3Meta(progressCh3, 'Field') = {field};
                        ch3Meta(progressCh3, 'ThisTilePath') = {ThisTilePath};
                    case 4
                        progressCh4 = progressCh4 + 1;
                        ch4Meta(progressCh4, 'StartRowThisTile') = {StartRowThisTile};
                        ch4Meta(progressCh4, 'StartColumnThisTile') = {StartColumnThisTile};
                        ch4Meta(progressCh4, 'Plane') = {plane};
                        ch4Meta(progressCh4, 'Field') = {field};
                        ch4Meta(progressCh4, 'ThisTilePath') = {ThisTilePath};
                end

            end
        end
    end
    ch1Ims = {};
    ch2Ims = {};
    ch3Ims = {};
    ch4Ims = {};

    parfor chi = 1:(height(ch1Meta))
        ch1Ims{chi} = imread(ch1Meta.ThisTilePath{chi});
        ch2Ims{chi} = imread(ch2Meta.ThisTilePath{chi});
        ch3Ims{chi} = imread(ch3Meta.ThisTilePath{chi});
        ch4Ims{chi} = imread(ch4Meta.ThisTilePath{chi});
    end
    
    ch1 = zeros(length(YVecThisOrganoid) * ImSize(1), length(XVecThisOrganoid) * ImSize(2), 25, 'uint16');
    for i = 1:height(ch1Meta)
        ch1(ch1Meta{i, 'StartRowThisTile'}:ch1Meta{i, 'StartRowThisTile'} + ImSize(1)-1, ch1Meta{i, 'StartColumnThisTile'}:ch1Meta{i, 'StartColumnThisTile'} + ImSize(2)-1, ch1Meta{i, 'Plane'}) = ch1Ims{i};
    end
    %vol(ch1)
    ch2 = zeros(length(YVecThisOrganoid) * ImSize(1), length(XVecThisOrganoid) * ImSize(2), 25, 'uint16');
    for i = 1:height(ch1Meta)
        ch2(ch2Meta{i, 'StartRowThisTile'}:ch2Meta{i, 'StartRowThisTile'} + ImSize(1)-1, ch2Meta{i, 'StartColumnThisTile'}:ch2Meta{i, 'StartColumnThisTile'} + ImSize(2)-1, ch2Meta{i, 'Plane'}) = ch2Ims{i};
    end
    %vol(ch2)
    ch3 = zeros(length(YVecThisOrganoid) * ImSize(1), length(XVecThisOrganoid) * ImSize(2), 25, 'uint16');
    for i = 1:height(ch1Meta)
        ch3(ch1Meta{i, 'StartRowThisTile'}:ch3Meta{i, 'StartRowThisTile'} + ImSize(1)-1, ch3Meta{i, 'StartColumnThisTile'}:ch3Meta{i, 'StartColumnThisTile'} + ImSize(2)-1, ch3Meta{i, 'Plane'}) = ch3Ims{i};
    end
    %vol(ch3)
    ch4 = zeros(length(YVecThisOrganoid) * ImSize(1), length(XVecThisOrganoid) * ImSize(2), 25, 'uint16');
    for i = 1:height(ch1Meta)
        ch4(ch1Meta{i, 'StartRowThisTile'}:ch4Meta{i, 'StartRowThisTile'} + ImSize(1)-1, ch4Meta{i, 'StartColumnThisTile'}:ch4Meta{i, 'StartColumnThisTile'} + ImSize(2)-1, ch4Meta{i, 'Plane'}) = ch4Ims{i};
    end
    %vol(ch4)
    OrganoidLoadTime = toc;
    disp(['Needed ', num2str(OrganoidLoadTime), ' seconds to load organoid images'])
    
    
    clear 'ch1Ims' 'ch2Ims' 'ch3Ims' 'ch4Ims' 
    chFindBestPlanes = uint32(ch1)+ uint32(ch2) + uint32(ch3) + uint32(ch4);
    chFindBestPlanesSummary = squeeze(sum(chFindBestPlanes,[1,2]));
    %figure; plot(chFindBestPlanesSummary)
    %GoodPlanes = chFindBestPlanesSummary > (max(chFindBestPlanesSummary)/5); % Above 20% of maximum
    GoodPlanes = chFindBestPlanesSummary > (max(chFindBestPlanesSummary)/3); % Above 33% of maximum
    GoodPlanes = bwareafilt(GoodPlanes, 1);
    clear 'chFindBestPlanes'
    ch1 = ch1(:,:,GoodPlanes);
    ch2 = ch2(:,:,GoodPlanes);
    ch3 = ch3(:,:,GoodPlanes);
    ch4 = ch4(:,:,GoodPlanes);
    ObjectsThisOrganoid = f_imageAnalysis(organoid, ch1, ch2, ch3, ch4, PreviewPath,Sample,ThisWellString);
    ObjectsAll{organoid} = ObjectsThisOrganoid;
    dataThisOrganoid = ObjectsThisOrganoid;
    save([OutPath, filesep, 'Slide_', ThisWellString, '_data_', num2str(organoid), '.mat'], 'dataThisOrganoid')
    writetable(dataThisOrganoid, [OutPath, filesep, 'Slide_', ThisWellString, '_data_', num2str(organoid), '.csv'], 'WriteRowNames', true)
end

Channelstable = MetaData.Channels{:}(:,(1:11));
LightSourcetable= MetaData.LightSources{:};
writetable(Channelstable, [OutPath, filesep, 'Channelstable.csv'], 'WriteRowNames', true) 
writetable(LightSourcetable, [OutPath, filesep, 'LightSourcetable.csv'], 'WriteRowNames', true)
data = vertcat(ObjectsAll{:});
save([OutPath, filesep, 'data_all.mat'], 'data')
writetable(data, [OutPath, filesep, 'data_all.csv'], 'WriteRowNames', true)


