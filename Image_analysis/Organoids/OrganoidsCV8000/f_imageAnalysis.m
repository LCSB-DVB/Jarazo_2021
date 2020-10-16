%function [outputArg1,outputArg2] = f_imageAnalysis(inputArg1,inputArg2)
function [ObjectsThisOrganoid] = f_imageAnalysis(Label, ch1, ch2, ch3, ch4, PreviewPath,Sample,ThisWellString)
%
    % the higher the 3rd number the less intense (more grey tones)
    % vol(ch1, 0, 15000) % HOECHST 33342
    % vol(ch2, 0, 25000) % Alexa 488 >>> Tuj1
    % vol(ch3, 0, 15000) % Alexa 568 >>> TH
    % vol(ch4, 0, 15000) % Alexa 647 >>> MAP2
    %% Initialize variables
    %keyboard
    
     %% Tuj1 (ch2)
    Tuj1_FT = zeros(size(ch2), 'double');

    parfor p=1:size(ch2, 3)
        Tuj1_FT(:,:,p) = f_LPF_by_FFT(ch2(:,:,p), 'Butterworth', [7,1], 0);
    end
    %vol(Tuj1_FT*100, 0, 1) % imtool(Tuj1_FT(:,:,5))
    Tuj1Mask = Tuj1_FT > 0.0025;% vol(Tuj1Mask, 0, 1)
    clear 'Tuj1_FT'
    Tuj1Mask = bwareaopen(Tuj1Mask, 1000);
    Tuj1Mask = medfilt3(Tuj1Mask);

    %vol(Tuj1Mask)  
    %% Organoid Mask (to remove the autofluorescence background)
    OrganoidMaskArea = imdilate(imdilate(Tuj1Mask, strel('disk', 50)), strel('sphere',5));
    
    %% Segment nuclei
    %vol(ch3, 0, 6000)
    ch1DoG = imfilter(ch1, fspecial('gaussian', 101, 3) - fspecial('gaussian', 101, 7), 'symmetric');%vol(ch1DoG, 0, 10)
    NucleiMask = ch1DoG > 30; %vol(NucleiMask)
    clear 'ch1DoG'
    NucleiMask = bwareaopen(NucleiMask, 300);%vol(NucleiMask)
    ch1LP = imfilter(ch1, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch1LP, 0, 15000, 'hot')
    NucMaskHigh =  (ch1LP > 4500) .* NucleiMask; %vol(NucMaskHigh, 0, 1)
    NucMaskAlive = NucleiMask & ~NucMaskHigh & OrganoidMaskArea; % vol(NucMaskAlive)
    
     %% TH (ch3)
    TH_FT = zeros(size(ch3), 'double');

    parfor p=1:size(ch3, 3)
        TH_FT(:,:,p) = f_LPF_by_FFT(ch3(:,:,p), 'Butterworth', [7,1], 0);
    end
    %vol(TH_FT*100, 0,1, 'hot')
    THMask = TH_FT > 0.0015;
    clear 'TH_FT'
    THMask = bwareaopen(THMask, 1000);
    THMask = THMask & OrganoidMaskArea;
    THMask =medfilt3(THMask);
    %vol(THMask)   
    %vol(TH_FT)    

    
    %% MAP2 (ch4)
    MAP2_FT = zeros(size(ch4), 'double');

    parfor p=1:size(ch4, 3)
        MAP2_FT(:,:,p) = f_LPF_by_FFT(ch4(:,:,p), 'Butterworth', [7,1], 0);
    end
    %vol(MAP2_FT*100, 0, 1) % imtool(MAP2_FT(:,:,5))
    MAP2Mask = MAP2_FT > 0.0025; %vol(MAP2Mask, 0, 1)
    clear 'MAP2_FT'
    MAP2Mask = bwareaopen(MAP2Mask, 1000);
    MAP2Mask = MAP2Mask & OrganoidMaskArea;
    MAP2Mask =medfilt3(MAP2Mask);    
   
    
%% TH skeleton3D EPFL
    disp('Start skel')
    tic
    skelTH = Skeleton3D(THMask);
    toc
    disp('Skel done')
%     vol(skelTH, 0, 1)
    [~, nodeTH, linkTH] = Skel2Graph3D(skelTH,0);                       

    ZeroNodeExplanationNeeded = 0;
    if ZeroNodeExplanationNeeded
        ZeroNodes = find(NodeDegreeVectorTH == 0);
        ZeroNodesLinIdx = vertcat(nodeTH(ZeroNodes).idx);
        ZeroNodeMask = zeros(size(THMaskClipped), 'uint8');
        ZeroNodeMask(ZeroNodesLinIdx) = 1; %vol(ZeroNodeMask)
        NodePreviewZeroCase = uint8(skelTH) + NodeMaskTH + 10*uint8(ZeroNodeMask) + uint8(THMask);
    end  
    
    %% TH Fragmentation
    
    % Define structuring element for surface detection
    Conn6 = strel('sphere', 1); % 6 connectivity
    % Detect surface
    SurfaceTH = THMask & ~(imerode(THMask, Conn6));
    %vol(SurfaceTH)  
    
    %% Perinuclear Volume 
    
    %vol(NucleiMask)
    NucleiMaskSingleCells = f_RemoveBigObjects (NucleiMask, 10000); 
    NucDil = imdilate(imdilate(NucleiMaskSingleCells, strel('disk', 4)), strel('sphere',1));
    NucPerim = logical(NucDil) & ~logical(NucleiMaskSingleCells);
    %vol(NucPerim)
   
    %% Percent TH pos
    %split perinuc
    D = bwdist(NucleiMaskSingleCells);
    %vol(D, 0, 20, 'hot')
    %it(D(:,:,1))
    disp('start watershed')
    tic
    W = watershed(D);
    toc
    disp('watershed done')
    %vol(W)
    NucPerimStencil = uint16(W) .* uint16(imreconstruct(logical(imdilate(NucPerim, strel('disk', 1))), logical(NucleiMaskSingleCells))); % This line was causing the error 20171207 % Function imreconstruct expected MARKER and MASK to have the same class.
    %vol(NucPerimStencil)
    %vol(NucPerim)
    %vol(NucleiMaskSingleCells)
       
    PeriNucMask = logical(NucPerimStencil);
    PeriNucMask = bwareaopen(PeriNucMask, 500);
    %vol(PeriNucMask)
    
    PerinucLM = bwlabeln(PeriNucMask);%vol(PerinucLM); vol(uint16(PeriNucMask) .* uint16(THMask), 0,1); vol(THMask +2*PeriNucMask)
    PeriNucObjects = regionprops('table', PerinucLM, double(THMask), 'PixelValues');
    THproportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjects, 'InputVariables', 'PixelValues');
    THPos = array2table(table2array(THproportions) > 0.01);
    THPos.Properties.VariableNames(end) = {'THpos'};
    PeriNucObjects = [PeriNucObjects, THproportions, THPos];
    PeriNucObjects.Properties.VariableNames(end-1) = {'THproportion'};%{'PixelValues', 'THproportion', 'THpos'};
    PeriNucObjectsCompact = PeriNucObjects(:, {'THproportion','THpos'});
    THPercent = (sum(PeriNucObjectsCompact.THpos)/height(PeriNucObjectsCompact))*100;
     
    PeriNucObjectsTuj1 = regionprops('table', PerinucLM, double(Tuj1Mask), 'PixelValues');
    Tuj1proportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjectsTuj1, 'InputVariables', 'PixelValues');
    Tuj1Pos = array2table(table2array(Tuj1proportions) > 0.01);
    Tuj1Pos.Properties.VariableNames(end) = {'Tuj1pos'};
    PeriNucObjectsTuj1 = [PeriNucObjectsTuj1, Tuj1proportions, Tuj1Pos];
    PeriNucObjectsTuj1.Properties.VariableNames(end-1) = {'Tuj1proportion'};
    PeriNucObjectsCompactTuj1 = PeriNucObjectsTuj1(:, {'Tuj1proportion','Tuj1pos'});
    Tuj1Percent = (sum(PeriNucObjectsCompactTuj1.Tuj1pos)/height(PeriNucObjectsCompactTuj1))*100;
    %% Previews 
    % Middle plane
    Midplane = round(size(ch1, 3)/2);
    % Scalebar
    imSize = [size(ch2, 1), size(ch2, 2)];
    [BarMask, ~] = f_barMask(200, 0.32393102760889064, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)

    PreviewTH = imadjust(ch3(:,:,Midplane));
    %imtool(PreviewTH)
    PreviewTHArea = imoverlay2(imadjust(ch3(:,:,Midplane)),THMask(:,:,Midplane), [1 0 0]);
    PreviewTHArea = imoverlay2(PreviewTHArea, BarMask, [1 1 1]);
    PreviewTHPerim = imoverlay2(imadjust(ch3(:,:,Midplane)), bwperim(THMask(:,:,Midplane)), [1 0 0]);
    PreviewTHPerim = imoverlay2(PreviewTHPerim, BarMask, [1 1 1]);
    
    PreviewNucAll = imadjust(ch1(:,:,Midplane));
    %imtool(PreviewNucAll)
    PreviewNucAllArea = imoverlay2(imadjust(ch1(:,:,Midplane)),NucleiMask(:,:,Midplane), [1 0 0]);
    PreviewNucAllArea = imoverlay2(PreviewNucAllArea, BarMask, [1 1 1]);
    PreviewNucAllPerim = imoverlay2(imadjust(ch1(:,:,Midplane)), bwperim(NucleiMask(:,:,Midplane)), [1 0 0]);
    PreviewNucAllPerim = imoverlay2(PreviewNucAllPerim, BarMask, [1 1 1]);
    
   
    PreviewNucAliveArea = imoverlay2(imadjust(ch1(:,:,Midplane)),NucMaskAlive(:,:,Midplane), [1 0 0]);
    PreviewNucAliveArea = imoverlay2(PreviewNucAliveArea, BarMask, [1 1 1]);
    PreviewNucAlivePerim = imoverlay2(imadjust(ch1(:,:,Midplane)), bwperim(NucMaskAlive(:,:,Midplane)), [1 0 0]);
    PreviewNucAlivePerim = imoverlay2(PreviewNucAlivePerim, BarMask, [1 1 1]);

    PreviewTuj1 = imadjust(ch2(:,:,Midplane));
    %imtool(PreviewTuj1)
    PreviewTuj1Area = imoverlay2(imadjust(ch2(:,:,Midplane)),Tuj1Mask(:,:,Midplane), [1 0 0]);
    PreviewTuj1Area = imoverlay2(PreviewTuj1Area, BarMask, [1 1 1]);
	PreviewTuj1Perim = imoverlay2(imadjust(ch2(:,:,Midplane)), bwperim(Tuj1Mask(:,:,Midplane)), [1 0 0]);
    PreviewTuj1Perim = imoverlay2(PreviewTuj1Perim, BarMask, [1 1 1]);
    
    PreviewMAP2 = imadjust(ch4(:,:,Midplane));
    %imtool(PreviewMAP2)
    PreviewMAP2Area = imoverlay2(imadjust(ch4(:,:,Midplane)),MAP2Mask(:,:,Midplane), [1 0 0]);
    PreviewMAP2Area = imoverlay2(PreviewMAP2Area, BarMask, [1 1 1]);
	PreviewMAP2Perim = imoverlay2(imadjust(ch4(:,:,Midplane)), bwperim(MAP2Mask(:,:,Midplane)), [1 0 0]);
    PreviewMAP2Perim = imoverlay2(PreviewMAP2Perim, BarMask, [1 1 1]);
    
    RGB_MAP2_1p = cat (3,  imadjust(ch3(:,:,Midplane)), imadjust(ch4(:,:,Midplane)), imadjust(ch1(:,:,Midplane)));
    RGB_MAP2_1p = imoverlay2(RGB_MAP2_1p, BarMask, [1 1 1]);
    RGB_MAP2_all = cat (3,  imadjust(max(ch3,[],3)), imadjust(max(ch4,[],3)), imadjust(max(ch1,[],3)));
    RGB_MAP2_all = imoverlay2(RGB_MAP2_all, BarMask, [1 1 1]);
    RGB_TUJ1_1p = cat (3,  imadjust(ch3(:,:,Midplane)), imadjust(ch2(:,:,Midplane)), imadjust(ch1(:,:,Midplane)));
    RGB_TUJ1_1p = imoverlay2(RGB_TUJ1_1p, BarMask, [1 1 1]);
    RGB_TUJ1_all = cat (3,  imadjust(max(ch3,[],3)), imadjust(max(ch2,[],3)), imadjust(max(ch1,[],3)));
    RGB_TUJ1_all = imoverlay2(RGB_TUJ1_all, BarMask, [1 1 1]);
    
    IdentityString = num2str(Label);
    imwrite(PreviewTH, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'TH_raw', '.png'])
    imwrite(PreviewTHArea, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'TH_Area', '.png'])
    imwrite(PreviewTHPerim, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'TH_Perim', '.png'])
    imwrite(PreviewNucAll, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'NucAll_raw', '.png'])
    imwrite(PreviewNucAllArea, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'NucAll_Area', '.png'])
    imwrite(PreviewNucAllPerim, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'NucAll_Perim', '.png'])
    imwrite(PreviewNucAliveArea, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'NucAlive_Area', '.png'])
    imwrite(PreviewNucAlivePerim, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'NucAlive_Perim', '.png'])
    imwrite(PreviewTuj1, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'Tuj1_raw', '.png'])
    imwrite(PreviewTuj1Area, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'Tuj1_Area', '.png'])
    imwrite(PreviewTuj1Perim, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'Tuj1_Perim', '.png'])
    imwrite(PreviewMAP2, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'MAP2_raw', '.png'])
    imwrite(PreviewMAP2Area, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'MAP2_Area', '.png'])
    imwrite(PreviewMAP2Perim, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'MAP2_Perim', '.png'])
    imwrite(RGB_MAP2_1p, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'RGB_MAP2_1p', '.png'])
    imwrite(RGB_MAP2_all, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'RGB_MAP2_all', '.png'])
    imwrite(RGB_TUJ1_1p, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'RGB_Tuj1_1p', '.png'])
    imwrite(RGB_TUJ1_all, [PreviewPath, filesep, Sample ,'_Slide_',ThisWellString,'_', IdentityString,'_' ,'RGB_Tuj1_all', '.png'])
   
    %% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = Label;
    ObjectsThisOrganoid.AreaName = convertCharsToStrings(Sample);
    ObjectsThisOrganoid.Slide = convertCharsToStrings(ThisWellString);
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.Tuj1MaskSum = sum(Tuj1Mask(:));
    ObjectsThisOrganoid.THMaskSum = sum(THMask(:));
    ObjectsThisOrganoid.MAP2MaskSum = sum(MAP2Mask(:)); 
    ObjectsThisOrganoid.MAP2ByTuj1 = sum(MAP2Mask(:)) / sum(Tuj1Mask(:));
    ObjectsThisOrganoid.Tuj1ByNuc = sum(Tuj1Mask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.Tuj1ByNucAlive = sum(Tuj1Mask(:)) / sum(NucMaskAlive(:));
    ObjectsThisOrganoid.MAP2ByNuc = sum(MAP2Mask(:))/ sum(NucleiMask(:));
    ObjectsThisOrganoid.MAP2ByNucAlive = sum(MAP2Mask(:))/ sum(NucMaskAlive(:));
    ObjectsThisOrganoid.THByNuc = sum(THMask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.THByNucAlive = sum(THMask(:)) / sum(NucMaskAlive(:));
    ObjectsThisOrganoid.THByTuj1 = sum(THMask(:)) / sum(Tuj1Mask(:));
    ObjectsThisOrganoid.THByMAP2 = sum(THMask(:)) / sum(MAP2Mask(:));
    ObjectsThisOrganoid.NucMaskDead = sum(NucMaskHigh(:));
    ObjectsThisOrganoid.NucMaskAlive = sum(NucMaskAlive(:));
    ObjectsThisOrganoid.THFragmentation = sum(SurfaceTH(:)) / sum(THMask(:));
    ObjectsThisOrganoid.SkelTH = sum(skelTH(:));
    ObjectsThisOrganoid.Nodes = size(nodeTH, 2);
    ObjectsThisOrganoid.Links = size(linkTH, 2);
    ObjectsThisOrganoid.THPercentofAllNuc = THPercent;
    ObjectsThisOrganoid.Tuj1PercentofAllNuc = Tuj1Percent;

end

