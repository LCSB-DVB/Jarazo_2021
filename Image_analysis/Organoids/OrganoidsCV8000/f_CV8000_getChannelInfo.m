function [MetaData] = f_CV8000_getChannelInfo(DataPath, MesPath)
% Load meta data from Yokogawa CV8000 MeasurementData folder
% Author: Paul Antony 20190607
    % Example:
    % DataPath = 'S:\YokogawaCV8000Horst\BTSData\CorrectedMeasurementData\Paul\PA neurons Training 20190326_20190326_153436\AssayPlate_Greiner_#781091';
    % MetaData = f_CV8000_getChannelInfo(DataPath)
    MetaData = table();

    %% Start from folder containing images
    MetaData.DataPath = DataPath;
    ImageFiles = dir([DataPath,filesep,'*.tif'])
    ImageFiles = struct2table(ImageFiles);
    
    %ImageFiles = [ImageFiles.folder, filesep, ImageFiles.name]
    ImageFiles = rowfun(@(a,b) {[a{:}, filesep, b{:}]}, ImageFiles, 'InputVariables', {'folder','name'})
    ImageFiles = table2cell(ImageFiles);
    %ImageFiles = sortrows(dirrec(DataPath, '.tif')');
    InfoTable = table();
    InfoTable.file = ImageFiles;
    Wells = regexp(InfoTable.file, '.*_(.\d{2})_T\d{4}F\d{3}L\d{2}.\d{2}Z\d{2}C\d{2}.tif', 'tokens');
    %% Remove calibration images
    SampleImBool = cellfun(@(x) ~isempty(x), Wells);
    Wells = Wells(SampleImBool);
    InfoTable = InfoTable(SampleImBool,:)
    InfoTable.Well = cellfun(@(x) x{:}{:}, Wells, 'UniformOutput', false);
    TimePoints = regexp(InfoTable.file, '.*_.\d{2}_T(\d{4})F\d{3}L\d{2}.\d{2}Z\d{2}C\d{2}.tif', 'tokens');
    InfoTable.Timepoint = cellfun(@(x) x{:}{:}, TimePoints, 'UniformOutput', false);
    Fields = regexp(InfoTable.file, '.*T\d{4}F(\d{3})L\d{2}.\d{2}Z\d{2}C\d{2}.tif', 'tokens');
    InfoTable.Field = cellfun(@(x) x{:}{:}, Fields, 'UniformOutput', false);
    TimeLines = regexp(InfoTable.file, '.*T\d{4}F\d{3}L(\d{2}).\d{2}Z\d{2}C\d{2}.tif', 'tokens');
    InfoTable.TimeLine = cellfun(@(x) x{:}{:}, TimeLines, 'UniformOutput', false);
    Planes = regexp(InfoTable.file, '.*T\d{4}F\d{3}L\d{2}.\d{2}Z(\d{2})C\d{2}.tif', 'tokens');
    InfoTable.Plane = cellfun(@(x) x{:}{:}, Planes, 'UniformOutput', false);
    Channels = regexp(InfoTable.file, '.*T\d{4}F\d{3}L\d{2}.\d{2}Z\d{2}C(\d{2}).tif', 'tokens');
    InfoTable.Channel = cellfun(@(x) x{:}{:}, Channels, 'UniformOutput', false);
    
    MetaData.InfoTable = {InfoTable};
    
    
    %% locate .mrf file and read bts:MeasurementSettingFileName, bts:OperatorName, bts:Title, bts:FieldCount, bts:ZCount
    MrfPath = dir([DataPath, filesep, '*.mrf']);
    MrfPath = [MrfPath.folder, filesep, MrfPath.name];
    %MrfPath = dirrec(DataPath, '.mrf'); MrfPath = MrfPath{:};
    MetaData.MrfPath = MrfPath;
    MrfXml = xmlread(MrfPath);
    MeasurementDetail =  MrfXml.getChildNodes.item(0);% disp(MeasurementDetail.getNodeName)
    
    MetaData.MeasurementSettingFileName = char(MeasurementDetail.getAttribute('bts:MeasurementSettingFileName'));% located in S:\YokogawaCV8000Horst\BTSData\MeasurementSetting\
    MetaData.OperatorName = char(MeasurementDetail.getAttribute('bts:OperatorName'));
    MetaData.Title = char(MeasurementDetail.getAttribute('bts:Title'));
    MetaData.FieldCount = char(MeasurementDetail.getAttribute('bts:FieldCount'));
    MetaData.ZCount = char(MeasurementDetail.getAttribute('bts:ZCount'));

    %% get .mes file, the measurement setting Xml

    %MesPath = ['S:\YokogawaCV8000Horst\BTSData\MeasurementSetting\', MetaData.MeasurementSettingFileName];
    %%MesPath = ['./', MetaData.MeasurementSettingFileName];
    %MesPath = '/mnt/irisgpfs/users/pantony/HCS_Projects/TemplateCV8000_20190703/IR_20190617_organoidRescanBinning1_25planes_0.mes'
    %%MetaData.MesPath = MesPath;
    
    MesXml = xmlread(MesPath);
    MeasurementSetting =  MesXml.getChildNodes.item(0);% disp(MeasurementSetting.getNodeName)

    
    

    %% get all time lines
    
    TimeLineList = MeasurementSetting.getElementsByTagName('bts:Timelapse').item(0);
    TimeLines = table();
    warning('off','MATLAB:table:RowsAddedExistingVars');
    TimeLineThis = TimeLineList.getFirstChild;
    TimeLineProgress = 0;
    while ~isempty(TimeLineThis)
        if ~max(TimeLineThis.getClass == 'class org.apache.xerces.dom.DeferredTextImpl')% skip #text
            TimeLineProgress = TimeLineProgress + 1;
            TimeLines(TimeLineProgress,'Name') = {char(TimeLineThis.getAttribute('bts:Name'))};
            TimeLines(TimeLineProgress,'InitialTime') = {char(TimeLineThis.getAttribute('bts:InitialTime'))};
            TimeLines(TimeLineProgress,'Period') = {char(TimeLineThis.getAttribute('bts:Period'))};
            TimeLines(TimeLineProgress,'Interval') = {char(TimeLineThis.getAttribute('bts:Interval'))};
            TimeLines(TimeLineProgress,'ExpectedTime') = {char(TimeLineThis.getAttribute('bts:ExpectedTime'))};
            %% Get well information
            WellsThisTimeline = TimeLineThis.getElementsByTagName('bts:WellSequence').item(0);
            TargetWells = table();
            TargetWell = WellsThisTimeline.getElementsByTagName('bts:TargetWell').item(0);
            TimelineWellProgress = 0;
            while ~isempty(TargetWell)
                if ~max(TargetWell.getClass == 'class org.apache.xerces.dom.DeferredTextImpl')% skip #text
                    TimelineWellProgress = TimelineWellProgress + 1;
                    TargetWells(TimelineWellProgress, 'WellRow') = {str2double(char(TargetWell.getAttribute('bts:Row')))};
                    TargetWells(TimelineWellProgress, 'WellCol') = {str2double(char(TargetWell.getAttribute('bts:Column')))};
                end
                TargetWell = TargetWell.getNextSibling;
            end
            TimeLines(TimeLineProgress,'Wells') = {TargetWells};
            %% Get field information
            FieldsThisTimeline = TimeLineThis.getElementsByTagName('bts:PointSequence').item(0);
            TargetFields = table();
            FieldSequence = FieldsThisTimeline.getElementsByTagName('bts:FixedPosition').item(0);
            FieldThis = FieldSequence.getElementsByTagName('bts:Point').item(0);
            TimelineFieldProgress = 0;
            while ~isempty(FieldThis)
                if ~max(FieldThis.getClass == 'class org.apache.xerces.dom.DeferredTextImpl')% skip #text
                    TimelineFieldProgress = TimelineFieldProgress + 1;
                    TargetFields(TimelineFieldProgress, 'Field') = {TimelineFieldProgress};
                    TargetFields(TimelineFieldProgress, 'X') = {str2double(char(FieldThis.getAttribute('bts:X')))};
                    TargetFields(TimelineFieldProgress, 'Y') = {str2double(char(FieldThis.getAttribute('bts:Y')))};
                end
                FieldThis = FieldThis.getNextSibling;
            end
            TimeLines(TimeLineProgress,'Fields') = {TargetFields};
            %%
        end

        TimeLineThis = TimeLineThis.getNextSibling;
    end
    MetaData.TimeLines = {TimeLines};
    
    
%     %% get all Wells
%     
%     Wellsequence = TimeLineList.getElementsByTagName('bts:Timelapse').item(0);
%     TimeLines = table();
%     warning('off','MATLAB:table:RowsAddedExistingVars');
%     TimeLineThis = TimeLineList.getFirstChild;
%     TimeLineProgress = 0;
%     while ~isempty(TimeLineThis)
%         if ~max(TimeLineThis.getClass == 'class org.apache.xerces.dom.DeferredTextImpl')% skip #text
%             TimeLineProgress = TimeLineProgress + 1;
%             TimeLines(TimeLineProgress,'Name') = {char(TimeLineThis.getAttribute('bts:Name'))};
%             TimeLines(TimeLineProgress,'InitialTime') = {char(TimeLineThis.getAttribute('bts:InitialTime'))};
%             TimeLines(TimeLineProgress,'Period') = {char(TimeLineThis.getAttribute('bts:Period'))};
%             TimeLines(TimeLineProgress,'Interval') = {char(TimeLineThis.getAttribute('bts:Interval'))};
%             TimeLines(TimeLineProgress,'ExpectedTime') = {char(TimeLineThis.getAttribute('bts:ExpectedTime'))};
%         end
%         TimeLineThis = TimeLineThis.getNextSibling;
%     end
%     MetaData.TimeLines = {TimeLines};
    
    %% get all light sources
    
    LightSourceList = MeasurementSetting.getElementsByTagName('bts:LightSourceList').item(0);
    LightSources = table();
    warning('off','MATLAB:table:RowsAddedExistingVars');
    LightSourceThis = LightSourceList.getFirstChild;
    LightSourceProgress = 0;
    while ~isempty(LightSourceThis)
        if ~max(LightSourceThis.getClass == 'class org.apache.xerces.dom.DeferredTextImpl')% skip #text
            LightSourceProgress = LightSourceProgress + 1;
            LightSources(LightSourceProgress,'LightSource') = {char(LightSourceThis.getAttribute('bts:Name'))};
            LightSources(LightSourceProgress, 'Power') = {char(LightSourceThis.getAttribute('bts:Power'))};
        end
        LightSourceThis = LightSourceThis.getNextSibling;
    end
    MetaData.LightSources = {LightSources};
    
    %% get all channels
    
    ChannelList = MeasurementSetting.getElementsByTagName('bts:ChannelList').item(0);
    Channels = table();
    warning('off','MATLAB:table:RowsAddedExistingVars');
    ChannelThis = ChannelList.getFirstChild;
    ChannelProgress = 0;
    while ~isempty(ChannelThis)
        if ~max(ChannelThis.getClass == 'class org.apache.xerces.dom.DeferredTextImpl')% skip #text
            ChannelProgress = ChannelProgress + 1;
            Channels(ChannelProgress,'ID') = {char(ChannelThis.getAttribute('bts:Ch'))};
            Channels(ChannelProgress, 'Objective') = {char(ChannelThis.getAttribute('bts:Objective'))};
            Channels(ChannelProgress, 'Objective') = {char(ChannelThis.getAttribute('bts:Objective'))};
            Channels(ChannelProgress, 'MethodID') = {char(ChannelThis.getAttribute('bts:Magnification'))};
            Channels(ChannelProgress, 'Method') = {char(ChannelThis.getAttribute('bts:MethodID'))};
            Channels(ChannelProgress, 'Method') = {char(ChannelThis.getAttribute('bts:Method'))};
            Channels(ChannelProgress, 'FilterID') = {char(ChannelThis.getAttribute('bts:FilterID'))};
            Channels(ChannelProgress, 'Acquisition') = {char(ChannelThis.getAttribute('bts:Acquisition'))};
            Channels(ChannelProgress, 'ExposureTime') = {char(ChannelThis.getAttribute('bts:ExposureTime'))};
            Channels(ChannelProgress, 'Binning') = {char(ChannelThis.getAttribute('bts:Binning'))};
            Channels(ChannelProgress, 'PinholeDiameter') = {char(ChannelThis.getAttribute('bts:PinholeDiameter'))};
            Channels(ChannelProgress, 'Kind') = {char(ChannelThis.getAttribute('bts:Kind'))};
            Channels(ChannelProgress, 'Fluorophore') = {char(ChannelThis.getAttribute('bts:Fluorophore'))};
            %Channels(ChannelProgress, 'Excitation') = {char(ChannelThis.getElementsByTagName('bts:LightSourceName').item(0).getTextContent())};
            %ExcitationThis = ChannelThis.getElementsByTagName('bts:LightSourceName').item(0);
            LightSourcesThisChannel = {};
            %warning('off','MATLAB:table:RowsAddedExistingVars');
            %LightSourceThis = ChannelThis.getFirstChild;
            LightSourceThis = ChannelThis.getElementsByTagName('bts:LightSourceName').item(0);
            LSProgress = 0;
            while ~isempty(LightSourceThis)
                if ~max(LightSourceThis.getClass == 'class org.apache.xerces.dom.DeferredTextImpl')% skip #text
                    LSProgress = LSProgress + 1;
                    LightSourcesThisChannel{LSProgress,1} = char(LightSourceThis.getTextContent);
                end
                LightSourceThis = LightSourceThis.getNextSibling;
            end
            Channels(ChannelProgress, 'Excitation') = {LightSourcesThisChannel};
        end
        ChannelThis = ChannelThis.getNextSibling;
    end
    MetaData.Channels = {Channels};
    
    %% Get coordinates from .icr
    %Coordinates = table();
    
    %IcrPath = dirrec(DataPath, '.icr');
    IcrPath = dir([DataPath, filesep, '*.icr']);
    %ImageFiles = dir([DataPath,filesep,'*.tif'])
    
    if size(IcrPath, 1) > 0 
        IcrPath = IcrPath{:};
        MetaData.IcrPath = IcrPath;
        IcrXml = xmlread(IcrPath);
        CoordinatesDetail =  IcrXml.getChildNodes.item(0);% disp(CoordinatesDetail.getNodeName)
    %     MeasurementBlockList = CoordinatesDetail.getElementsByTagName('ict:MeasurementBlockList').item(0);
    %     MeasurementBlockList.getChildNodes.item(1)
        WellResultList = CoordinatesDetail.getElementsByTagName('ict:WellResultList').item(0);
        %WellCount = WellResultList.getLength
        %% Loop over wells
        %%%%%%%%%%%%%%%%%
        Coordinates = table();
        WellThis = WellResultList.getFirstChild;
        ImageProgress = 0;
        while ~isempty(WellThis)
            if ~max(WellThis.getClass == 'class org.apache.xerces.dom.DeferredTextImpl')% skip #text
                ImageResultThis = WellThis.getElementsByTagName('ict:ImageResult').item(0);% Need to loop
                while ~isempty(ImageResultThis)
                    if ~max(ImageResultThis.getClass == 'class org.apache.xerces.dom.DeferredTextImpl')% skip #text
                        ImageProgress = ImageProgress + 1;
                        Coordinates(ImageProgress,'MeasurementBlockIndex') = {str2double(char(ImageResultThis.getAttribute('ict:MeasurementBlockIndex')))};
                        %Coordinates(ImageProgress,'ROW') = {char(ImageResultThis.getAttribute('ict:Row'))};
                        Coordinates(ImageProgress,'ROW') = {str2double(char(ImageResultThis.getAttribute('ict:Row')))};
                        Coordinates(ImageProgress,'COLUMN') = {str2double(char(ImageResultThis.getAttribute('ict:Column')))};
                        Coordinates(ImageProgress,'FIELD') = {str2double(char(ImageResultThis.getAttribute('ict:Field')))};
                        Coordinates(ImageProgress,'ZIndex') = {str2double(char(ImageResultThis.getAttribute('ict:ZIndex')))};
                        Coordinates(ImageProgress,'TimePoint') = {str2double(char(ImageResultThis.getAttribute('ict:TimePoint')))};
                        Coordinates(ImageProgress,'TimeLine') = {str2double(char(ImageResultThis.getAttribute('ict:TimeLine')))};
                        Coordinates(ImageProgress,'ActionIndex') = {str2double(char(ImageResultThis.getAttribute('ict:ActionIndex')))};
                        Coordinates(ImageProgress,'Channel') = {str2double(char(ImageResultThis.getAttribute('ict:Channel')))};
                        Coordinates(ImageProgress,'ShiftX') = {str2double(char(ImageResultThis.getAttribute('ict:ShiftX')))};
                        Coordinates(ImageProgress,'ShiftY') = {str2double(char(ImageResultThis.getAttribute('ict:ShiftY')))};
                        Coordinates(ImageProgress,'IsSaturationPixelFromSourceImage') = {char(ImageResultThis.getAttribute('ict:IsSaturationPixelFromSourceImage'))};
                        Coordinates(ImageProgress,'OverLimitMaxIntensity') = {char(ImageResultThis.getAttribute('ict:OverLimitMaxIntensity'))};
                    end
                    ImageResultThis = ImageResultThis.getNextSibling;
                end
            end
            WellThis = WellThis.getNextSibling;
        end
        MetaData.Coordinates = {Coordinates};
        %%%%%%%%%%%%%%%%%
    % %     %WellResult = WellResultList.getElementsByTagName('ict:WellResult').item(0);% Need to loop
    % %     %WellResultList = MeasurementBlockList.getElementsByTagName('ict:WellResultList').item(0);
    % %     ImageResultList = WellResult.getElementsByTagName('ict:ImageResultList').item(0);% Need to loop
    % %     ImageResult = ImageResultList.getChildNodes.item(0)
    % %     Coordinates.MeasurementBlockIndex
    % %     Coordinates.ROW
    % %     Coordinates.COLUMN
    % %     Coordinates.FIELD
    % %     Coordinates.ZIndex
    % %     Coordinates.TimePoint
    % %     Coordinates.TimeLine
    % %     Coordinates.ActionIndex
    % %     Coordinates.Channel
    % %     Coordinates.ShiftX
    % %     Coordinates.ShiftY
    % %     Coordinates.IsSaturationPixelFromSourceImage
    % %     Coordinates.OverLimitMaxIntensity
    end
end