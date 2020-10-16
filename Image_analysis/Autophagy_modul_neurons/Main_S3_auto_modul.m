%% Clear Matlab workspace
clear
clc

%% User inputs per run

Barcodes = {'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200728_10X_Autophagy_Modulation_Chloro_S3_C2_R1',...
		'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200728_10X_Autophagy_Modulation_Chloro_S3_C2_R2',...
		'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200728_10X_Autophagy_Modulation_Chloro_S3_C2_R3'};


for b = 1:size(Barcodes, 2) %loop over barcodes
    Barcode = Barcodes{b};
    Objects = {};   
    clear ObjectsAll
    ObjectsCell = {};
    
    %% Document script
    BarcodeName = regexp(Barcode, '.*\\(.*)', 'tokens'); 
    BarcodeName = BarcodeName{:};
    SavePathBegining = {'Y:\16-Our Papers\In Preparation\PINK1\PNAS\Revision\Figures\Figure_3\Originals\I\Script_output\'};
    SavePathMainDirectory = [SavePathBegining{:}, BarcodeName{:}];
    AnalysisTimeStamp = datestr(now, 'yyyymmdd_HHMMSS');
    SavePath = [SavePathMainDirectory, '\Analysis_', AnalysisTimeStamp];
    mkdir(SavePath)
    FileNameShort = mfilename;
    newbackup = sprintf('%s_log.m',[SavePath, '\', FileNameShort]);
    FileNameAndLocation = mfilename('fullpath');
    currentfile = strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup); 
    PreviewPath = [SavePath, filesep, 'Previews'];
    mkdir(PreviewPath)
    Version = version();
    save([SavePath filesep 'MatlabVersion.mat'], 'Version')

    f_LogDependencies(FileNameShort, SavePath); 
    %% Load data
    InfoTable = f_InfoTable(Barcodes{b});

    parfor f=1:height(InfoTable)%loop over fields
        InfoTableThis = InfoTable(f,:);
        cube = readflexcube(InfoTable.files{f}, 'PlaneCount', 1); % Read 4-D image cube
        ch1 = cube.data(:,:,:,1); %Hoechst it(ch1)
        ch2 = cube.data(:,:,:,2); %Tuj1 % it(ch2)
        ch3 = cube.data(:,:,:,3); %TH % it(ch3)
        ch4 = cube.data(:,:,:,4); %GFAP % it(ch4)
       
        Objects = ImageAnalysis_S3_auto_modul(InfoTableThis,SavePath, ch1, ch2, ch3, ch4);
        
        if size(Objects,1) == 1
           ObjectsCell{f} = Objects;
        end
            
    end

      ObjectsAll = vertcat(ObjectsCell{:});
      save([SavePath, filesep, 'data.mat'], 'ObjectsAll');
      writetable(ObjectsAll, [SavePath, '\Objects.csv'], 'WriteRowNames', true)  
      writetable(ObjectsAll, [SavePath, '\Objects.xls'], 'WriteRowNames', true)
      writetable(ObjectsAll, [SavePath, '\Objects.txt'])    
end

