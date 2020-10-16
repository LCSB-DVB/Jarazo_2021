function [InfoTable, MissingSamples] = f_InfoTable(BarcodePath, varargin)
%This function returns a table withh image paths, areanames, and well coordinates
%   BarcodePath: Full path to barcode in the OperaDB on LCSB_HCS
%   By Paul Antony 20161130
%   In case the layout is not found, provide the full path of the layout as
%   second input
%   Example of use:  InfoTable = f_InfoTable(Barcodes{b}, 'S:\OperaQEHS\OperaDB\Jonathan Arias\JA20161203_AllBasal_AutoMito_LC3\Settings\basal_conditions.lay');


    %% Read .apr file

    Barcode = regexp(BarcodePath, '\\.*\\(.*)', 'tokens'); Barcode = Barcode{:}{:};
    Aprfiles = dirrec(BarcodePath, '.apr');
    AprFile = Aprfiles{1}; % take first best apr file assuming that the layout does not change between Meas folders

    fileID = fopen(AprFile);

    tline = fgetl(fileID);

    while ischar(tline) % loop line by line until end of file
        % this loop assumes that there is only one type of .exs (.exp and .lay are unique by machine design)

        TokenLay = regexp(tline, '.*(OperaDB.*\.lay).*', 'tokens');
        if size(TokenLay,1) > 0
            Lay = TokenLay{:}{:}; % path to layout xml
            break
        end

        tline = fgetl(fileID);

    end

    fclose(fileID);

    %% Extract information fron the .lay layout XML
    
    if nargin == 1
        theNode = xmlread(['S:\OperaQEHS\', Lay]);  % XML
    elseif nargin == 2
        theNode = xmlread(varargin{1});  % XML
    end
    
    theNode = theNode.getChildNodes.item(0);    % Plate layout
    theNode = theNode.getChildNodes.item(2);    % The 3rd child node of Plate layout contains RectangularElements
    MetaData = cell(theNode.getLength,3);
    well = 0; % Progress variable
    
    for i = 1:theNode.getLength

        theRectangle = theNode.getChildNodes.item(i-1); % Rectangles contain AreaName information
        theAreaName = char(theRectangle.getChildNodes.item(0).getTextContent);
        AttributesThisNode = theRectangle.getAttributes;
        bottom = str2double(char(AttributesThisNode.item(0).getTextContent));
        left = str2double(char(AttributesThisNode.item(1).getTextContent));
        right = str2double(char(AttributesThisNode.item(2).getTextContent));
        top = str2double(char(AttributesThisNode.item(3).getTextContent));

        for row = top:bottom
            for col = left:right
                well = well + 1;
                MetaData{well,1} = row;
                MetaData{well,2} = col;
                MetaData{well,3} = {theAreaName};
            end
        end

    end

    MetaData = cell2table(MetaData);
    MetaData.Properties.VariableNames = {'Row','Column','AreaName'};


    %% Get full paths to images via OS

    Inprefix = 'S:\OperaQEHS\OperaArchiveCol';
    inpath = [Inprefix filesep Barcode filesep];
    files = dirrec(inpath, '.flex')';
    RowColField = regexp(files, '.*\\(.*).flex', 'tokens');
    Row = cellfun(@(x) str2double(x{:}{:}(1:3)), RowColField, 'UniformOutput', false);
    Column = cellfun(@(x) str2double(x{:}{:}(4:6)), RowColField, 'UniformOutput', false);
    field = cellfun(@(x) str2double(x{:}{:}(7:9)), RowColField, 'UniformOutput', false);
    files = [cell2table(files), cell2table(Row), cell2table(Column), cell2table(field)];

    %% Join layout and OS tables
    try
        InfoTable = join(files, MetaData);
        MissingSamples = 'There are no missing samples'
    catch
        InfoTable = outerjoin(files, MetaData);
        MissingSamples = InfoTable(isnan(InfoTable.Row_files),:);
        InfoTable = InfoTable(~isnan(InfoTable.Row_files),:);
        InfoTable.Properties.VariableNames = strrep(InfoTable.Properties.VariableNames, '_files', '');
    end
    FullPaths = InfoTable.files;
    %regexp('.*\\(.*)\\Meas', FullPaths, 'tokens')
    regexp('S:\OperaQEHS\OperaArchiveCol\JA_20161127_WTcontrol_AutoSensor\Meas_01(2016-11-27_21-09-43)\002003001.flex', '(.*)' , 'tokens');
    
    BarCodeList = cellfun(@(x) regexp(x, '.*\\(.*)\\Meas', 'tokens'), FullPaths, 'UniformOutput', false);
    BarCodeList = cellfun(@(x) x{:}{:}, BarCodeList, 'UniformOutput', false);
    
    InfoTable = [InfoTable, BarCodeList];
    InfoTable.Properties.VariableNames{end} = 'Barcode';

    
end

