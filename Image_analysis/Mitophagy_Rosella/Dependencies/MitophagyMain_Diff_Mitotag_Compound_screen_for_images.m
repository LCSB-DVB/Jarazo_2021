clear
clc

%Barcodes = {'JA_20161127_WTcontrol_MitoSensor'};
%Barcodes = 'JJ_20170824_hiPSC_Mitotag_compound_titration_N2';%,...
Barcodes = {'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171206_60X_NESC_Mito_Lyso_tags_DIFF_1\JJ_20171206_60X_NESC_Mito_Lyso_tags_DIFF_1'...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171206_60X_NESC_Mito_Lyso_tags_DIFF_2\JJ_20171206_60X_NESC_Mito_Lyso_tags_DIFF_2',...
    'S:\OperaQEHS\OperaDB\Javier Jarazo\JJ_20171206_60X_NESC_Mito_Lyso_tags_DIFF_3\JJ_20171206_60X_NESC_Mito_Lyso_tags_DIFF_3'};

%% Image analysis
Objects_Loop = {};

for b = 1:size(Barcodes,2)
    Objects = {};   
    clear Objects_Loop
    
    % Automated documentation
    Barcode = Barcodes{b};
    BarcodeName = regexp(Barcode, '.*\\(.*)', 'tokens'); 
    BarcodeName = BarcodeName{:};
    SavePathBegining = {'S:\HCS_Platform\Data\JavierJarazo\Diff_Mitotag_Compscreen\ATP5C1\'};
    SavePathMainDirectory = [SavePathBegining{:}, BarcodeName{:}];
    AnalysisTimeStamp = datestr(now, 'yyyymmdd_HHMMSS');
    FolderThisAnalysis = [SavePathMainDirectory, '\Analysis_', AnalysisTimeStamp];
    mkdir(FolderThisAnalysis)
    PreviewPath = [FolderThisAnalysis, '\Previews'];
    mkdir(PreviewPath);
    FileNameShort = mfilename;
    newbackup = sprintf('%s_log.m',[FolderThisAnalysis, '\', FileNameShort]);
    FileNameAndLocation = [mfilename('fullpath')];
    currentfile = strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup);
    Version = version;
    save([FolderThisAnalysis, filesep, 'Version.mat'], 'Version')
    f_LogDependencies(FileNameShort, FolderThisAnalysis); 
    
    InfoTable = f_InfoTable_Time_series(Barcodes{b});
     %[files, info, areaNames, plateLayout, areaMap] = f_Local_Read_Files_Timeseries_Mosaic_BT2('S:\OperaQEHS\OperaArchiveCol', Barcode, 1);
%      InfoTable = sortrows(struct2table(info), {'row', 'column', 'field'});
    %parfor i = 1:height(InfoTable)
    for i = 1:height(InfoTable)
        
        [Objects] = Mitophagy_compound_screen_for_images(InfoTable.files, i, InfoTable, PreviewPath);
        Objects_Loop{b,i} = Objects;

    end
    
    Objects = vertcat(Objects_Loop{:}); 
    saveTableToHDF5([FolderThisAnalysis, '\Objectshd5.mat'], 'ObjectList', Objects);
    save([FolderThisAnalysis, '\Objects.mat'], 'Objects');
    saveTableToHDF5([FolderThisAnalysis, '\Objects.hd5'], 'ObjectList', Objects);
    writetable(Objects, [FolderThisAnalysis, '\Objects.txt'])  
end





%%