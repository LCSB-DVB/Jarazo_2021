%% Clear Matlab workspace
clear
clc

%% User inputs per run

Barcodes = {'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d7_N1\JJ_20171219_20X_Diff_Apop_d7_N1',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d7_N2\JJ_20171219_20X_Diff_Apop_d7_N2',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d7_N3\JJ_20171219_20X_Diff_Apop_d7_N3',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d14_N1\JJ_20171219_20X_Diff_Apop_d14_N1',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d14_N2\JJ_20171219_20X_Diff_Apop_d14_N2',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d14_N3\JJ_20171219_20X_Diff_Apop_d14_N3',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d21_N1\JJ_20171219_20X_Diff_Apop_d21_N1',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d21_N2\JJ_20171219_20X_Diff_Apop_d21_N2',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d21_N3\JJ_20171219_20X_Diff_Apop_d21_N3',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171219_20X_Diff_Apop_d42\JJ_20171219_20X_Diff_Apop_d42'};


for b = 1:size(Barcodes, 2) %loop over barcodes
    Barcode = Barcodes{b};
    Objects = {};   
    clear ObjectsAll
    ObjectsCell = {};
    
    %% Document script
    % SavePathMainDirectory = ['.\', Barcode];
    BarcodeName = regexp(Barcode, '.*\\(.*)', 'tokens'); 
    BarcodeName = BarcodeName{:};
    SavePathBegining = {'S:\HCS_Platform\Data\JavierJarazo\Staging\Diff_'};
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
        ch2 = cube.data(:,:,:,2); %PARP% it(ch2)
        ch3 = cube.data(:,:,:,3); %TH % it(ch3)
        ch4 = cube.data(:,:,:,4); %Tuj1 % it(ch4)
       
        Objects = Analysis_Diff(InfoTableThis, f, SavePath, ch1, ch2, ch3, ch4);
        
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

