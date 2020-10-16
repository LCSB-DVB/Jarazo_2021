%% Mitochondrial Morphometrics Opera 3D 
clear
clc

%% User inputs

Barcodes = {'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200717_40X_Autophagy_Modulation_Bafi_S1_C2_R1',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200717_40X_Autophagy_Modulation_Bafi_S1_C2_R2',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200717_40X_Autophagy_Modulation_Bafi_S1_C2_R3',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200717_40X_Autophagy_Modulation_Chloro_S1_C2_R1',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200717_40X_Autophagy_Modulation_Chloro_S1_C2_R2',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200717_40X_Autophagy_Modulation_Chloro_S1_C2_R3',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200717_40X_Autophagy_Modulation_Rapa_S1_C2_R1',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200717_40X_Autophagy_Modulation_Rapa_S1_C2_R2',...
	'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20200717_40X_Autophagy_Modulation_Rapa_S1_C2_R3'};


for b = 1:size(Barcodes, 2) %loop over barcodes
    Barcode = Barcodes{b};
    Objects = {};   
    clear ObjectsAll
    ObjectsCell = {};
    
    %% Document script
    BarcodeName = regexp(Barcode, '.*\\(.*)', 'tokens'); 
    BarcodeName = BarcodeName{:};
    SavePathBegining = {'S:\HCS_Platform\Data\JavierJarazo\Autophagy_modulation\S1\'};
    SavePathMainDirectory = [SavePathBegining{:}, BarcodeName{:}];
    AnalysisTimeStamp = datestr(now, 'yyyymmdd_HHMMSS');
    SavePath = [SavePathMainDirectory, '\Analysis_', AnalysisTimeStamp];
    PreviewPath = [SavePath, filesep, 'Previews'];
    mkdir(SavePath)
    mkdir(PreviewPath)
    FileNameShort = mfilename;
    newbackup = sprintf('%s_log.m',[SavePath, '\', FileNameShort]);
    FileNameAndLocation = mfilename('fullpath');
    currentfile = strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup); 
    Version = version();
    save([SavePath filesep 'MatlabVersion.mat'], 'Version')

    f_LogDependencies(FileNameShort, SavePath); 
    %% Load data
    InfoTable = f_InfoTable(Barcodes{b});
    parfor l = 1:height(InfoTable)
            cube = readflexcube(InfoTable.files{l}, 'PlaneCount', 1); % Read 4-D image cube
            ch1 = cube.data(:,:,:,1);%Hoechst
            ch2 = cube.data(:,:,:,2);%LC3
            ch3 = cube.data(:,:,:,3);%LAMP1
			ch4 = cube.data(:,:,:,4);%TH % it(ch4)
            %vol(ch1)
            %vol(ch2)
            %vol(ch3)
            InfoTableThis = InfoTable(l,:);
            Objects = ImageAnalysis_S1_auto_modul(ch1,ch2,ch3,ch4,InfoTableThis, PreviewPath);
        if size(Objects,1) == 1
           ObjectsCell{l} = Objects;
        end
            
    end

      ObjectsAll = vertcat(ObjectsCell{:});
      save([SavePath, filesep, 'data.mat'], 'ObjectsAll');
      writetable(ObjectsAll, [SavePath, '\Objects.csv'], 'WriteRowNames', true)  
      writetable(ObjectsAll, [SavePath, '\Objects.xls'], 'WriteRowNames', true)
      writetable(ObjectsAll, [SavePath, '\Objects.txt'])
   
end

