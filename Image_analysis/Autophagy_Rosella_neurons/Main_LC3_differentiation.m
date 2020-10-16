%% Prepare Matlab
clear
clc

%% User inputs
Barcodes = {'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S1_C1_R1\JJ_20200810_60X_Auto_S1_C1_R1'...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S1_C1_R2\JJ_20200810_60X_Auto_S1_C1_R2',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S1_C1_R3\JJ_20200810_60X_Auto_S1_C1_R3',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S1_C2_R1\JJ_20200810_60X_Auto_S1_C2_R1',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S1_C2_R2\JJ_20200810_60X_Auto_S1_C2_R2',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S1_C2_R3\JJ_20200810_60X_Auto_S1_C2_R3',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S2_C1_R1\JJ_20200810_60X_Auto_S2_C1_R1'...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S2_C1_R2\JJ_20200810_60X_Auto_S2_C1_R2',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S2_C1_R3\JJ_20200810_60X_Auto_S2_C1_R3',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S2_C2_R1\JJ_20200810_60X_Auto_S2_C2_R1',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S2_C2_R2\JJ_20200810_60X_Auto_S2_C2_R2',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200810_60X_Auto_S2_C2_R3\JJ_20200810_60X_Auto_S2_C2_R3'};

Objects_Loop = {};

%% Loop over Barcodes and images
for b = 1:size(Barcodes,2)
	 
    clear Objects_Loop
% Automated documentation
    Barcode = Barcodes{b};
    BarcodeName = regexp(Barcode, '.*\\(.*)', 'tokens'); 
    BarcodeName = BarcodeName{:};
    SavePathBegining = 'S:\HCS_Platform\Data\JavierJarazo\Diff_Autotag_Compscreen\LC3\Revision\';
    SavePathMainDirectory = [SavePathBegining, BarcodeName{:}];
    AnalysisTimeStamp = datestr(now, 'yyyymmdd_HHMMSS');
    FolderThisAnalysis = [SavePathMainDirectory, '\Analysis_', AnalysisTimeStamp];
    mkdir(FolderThisAnalysis)
    PreviewPath = [FolderThisAnalysis, '\Previews'];
    mkdir(PreviewPath);
    FileNameShort = mfilename;
    newbackup = sprintf('%s_log.m',[FolderThisAnalysis, '\', FileNameShort]);
    FileNameAndLocation = mfilename('fullpath');
    currentfile = strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup);
    Version = version;
    save([FolderThisAnalysis, filesep, 'Version.mat'], 'Version')
    f_LogDependencies(FileNameShort, FolderThisAnalysis); 
    
    InfoTable = f_InfoTable_Time_series_ROW(Barcodes{b});    
 
     ObjectsAllFields = {};
	
        parfor i = 1:height(InfoTable)
            cube = readflexcube(InfoTable.files{i}, 'PlaneCount', 1); % Read 4-D image cube
            ch1 = cube.data(:, :, :, 1); %pHluorin (Green) %vi(flip(ch1,1))
            ch2 = cube.data(:, :, :, 2); %DsRed %vi(flip(ch2,1))
            InfoTableThis = InfoTable(i,:);
            [ObjectsAllFields]= ImageAnalysisLC3_differentiation(ch1,ch2,InfoTableThis, PreviewPath);
			Objects_Loop{b,i} = ObjectsAllFields;
        end % images
	
	Objects = vertcat(Objects_Loop{:});
	save([FolderThisAnalysis, '\Objects.mat'], 'Objects');
	writetable(Objects, [FolderThisAnalysis, '\Objects.txt'])
	
	ObjectsGrouped = grpstats(Objects, {'Barcode','AreaName', 'Column', 'Row', 'Field','TimePoint','Areaofcells'}, {'sum','mean','sem','std','median'});
    writetable(ObjectsGrouped, [FolderThisAnalysis,'\ObjectsGrouped.csv'])
    save([FolderThisAnalysis,'\SummaryAll.mat'], 'ObjectsGrouped');
	
end % Barcodes


