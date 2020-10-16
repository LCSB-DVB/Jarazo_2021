%% Prepare Matlab

% let the launcher do...

%% User inputs
decision3D = 1;
SizeExclusionThreshold = 7;
showIllustration = 0;
Barcodes = {Barcode};

CreateNewPSF = 0;


%% Initializations
Conn6Strel = strel('sphere', 1); Conn6Strel = Conn6Strel.Neighborhood;   
progress = 0;
load('FlatFieldIm.mat') % vol(CalImGreen)

%% Self documentation

AnalysisTimeStamp = datestr(now, 'yyyymmdd_HHMMSS');
SavePath = ['S:\HCS_Platform\Data\JavierJarazo\Rosella\LC3\Analysis_', AnalysisTimeStamp];
mkdir(SavePath)
FileNameShort = mfilename;
newbackup = sprintf('%s_log.m',[SavePath, '\', FileNameShort]);
FileNameAndLocation = mfilename('fullpath');
currentfile = strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);
Version = version();
save([SavePath filesep 'MatlabVersion.mat'], 'Version')

Dependencies = f_LogDependencies([FileNameShort, '.m'], SavePath);


%% Prepare output structure

% Create preview dir

if exist([SavePath filesep 'Previews'], 'dir') ~= 7
    mkdir([SavePath filesep 'Previews'])
end


% %% Load data

%% Loop over Barcodes and images
for b = 1:size(Barcodes,2)
    if exist('Layout')
        InfoTable = f_InfoTable(Barcodes{b}, Layout);
    else
        InfoTable = f_InfoTable(Barcodes{b});
    end
    
    files = InfoTable.files;
    
    
    for i = 1:height(InfoTable)
        progress = progress + 1;

        cube = readflexcube(files{i}, 'PlaneCount', 5); % Read 4-D image cube
        ch1 = cube.data(:, :, :, 1); %pHluorin (Green)
        %imwrite(imadjust(max(ch1,[],3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\pHluorin_raw.png')
        %%imwrite(imadjust(ch1(:,:,3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\pHluorin_raw.png')
        %it(imadjust(max(ch1,[],3)))
        ch2 = cube.data(:, :, :, 2); %DsRed
        if ismember(Barcodes{b}, {'S:\OperaQEHS\OperaDB\Jonathan Arias\JA_20170225_autolines_001','S:\OperaQEHS\OperaDB\Jonathan Arias\JA_20170225_autolines_002','S:\OperaQEHS\OperaDB\Jonathan Arias\JA_20170225_autolines_003',...
                                  'S:\OperaQEHS\OperaDB\Jonathan Arias\JA_20170225_autolines_004','S:\OperaQEHS\OperaDB\Jonathan Arias\JA_20170225_autolines_005'})
            ch1 = ApplyFlatFieldCorrection(ch1, CalImGreen); % vol(ch1)
            ch2 = ApplyFlatFieldCorrection(ch2, CalImRed); % vol(ch2)
        end
        %imwrite(imadjust(max(ch2,[],3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\DsRed_raw.png')
        %%imwrite(imadjust(ch2(:,:,3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\DsRed_raw.png')
        %vi(flip(ch1,1)) % vol(ch1)
        %vi(flip(ch2,1))

        %% Deconvolution
        [ImR, ImG] = Deconvolution_Javier(CreateNewPSF, ch2, ch1);
        %imwrite(imadjust(max(ImR,[],3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\ImRdeconvolved.png')
        %imwrite(imadjust(max(ImG,[],3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\ImGdeconvolved.png')
        %%imwrite(imadjust(ImR(:,:,3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\ImRdeconvolved.png')
        %%imwrite(imadjust(ImG(:,:,3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\ImGdeconvolved.png')
        %vi(flip(ImG,1)) vol(ImG, 0, 10000)
        % toto = adapthisteq(ImG(:,:,3),'NumTiles', [500,500]);
        % vol(toto) % Do not uncoment, out of of focus is real
        %vi(flip(ImR,1))

        %% Segment all LC3 objects DsRed

        RedallDoG = imfilter(ImR, fspecial('gaussian', 20, 1), 'symmetric') - imfilter(ImR, fspecial('gaussian', 20, 7), 'symmetric'); 
        %imwrite(imadjust(max(RedallDoG,[],3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\RedAllDoG.png')
        %%imwrite(imadjust(RedallDoG(:,:,3)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\RedAllDoG.png')
        %vi(flip(RedallDoG,1))
        %vol(RedallDoG, 0, 5000, 'hot')
        RedallMask = RedallDoG > 400; %vi(flip(uint8(RedallMask),1))
       

        
        RedTophat = imtophat(ImR, strel('disk', 25)); %vol(RedTophat, 0, 10000, 'hot')
        RedTophatMask = RedTophat > 1200; %vol(RedTophatMask)
        %vol(RedTophatMask +2*RedallMask)
        SizeThreshold = 500;
        ProportionThreshold = 0.1;
        MaskReconstructed = SegmentFilterSmallInBig(RedTophatMask, RedallMask, SizeThreshold, ProportionThreshold); % vol(MaskReconstructed)
        RedallMask = (RedallMask + MaskReconstructed) > 0; % vol(RedallMask)
        
        
        RedallMask = bwareaopen(RedallMask, SizeExclusionThreshold); % vol(RedallMask)
        %imwrite(imadjust(uint16(max(RedallMask,[],3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\smallDsRedvesicles.png')
        %%imwrite(imadjust(uint16(RedallMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\smallDsRedvesicles.png')

        %RedBigDoG = imfilter(ImR, fspecial('gaussian', 60, 1), 'symmetric') - imfilter(ImR, fspecial('gaussian', 60, 20), 'symmetric');
        RedBigDoG = imfilter(ImR, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ImR, fspecial('gaussian', 25, 6), 'symmetric');
        %vi(flip(RedBigDoG,1)) % vol(RedBigDoG, 0,10000, 'hot')
        RedBigMask = RedBigDoG > 1000; %vi(flip(uint8(RedBigMask),1))
        RedBigMask = bwareaopen(RedBigMask, 200);
        RedBigMask = f_RemoveBigObjects(RedBigMask, 2000);
        RedBigMaskPadded = padarray(RedBigMask,[0 0 1]);
        RedBigMaskPadded = imclearborder(RedBigMaskPadded);
        RedBigMask = RedBigMaskPadded(:,:,2:end-1); % vol(RedBigMask)
        %imwrite(imadjust(uint16(max(RedBigMask,[],3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\bigDsRedvesicles.png')
        %%imwrite(imadjust(uint16(RedBigMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\bigDsRedvesicles.png')

        % split vesicles
        RedDist = bwdist(RedallMask);% Use small objects as seed
        %it(imadjust(uint16(max(RedDist,[],3))))
        %imwrite(imadjust(uint16(max(RedDist,[],3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\RedDist.png')
        %%imwrite(imadjust(uint16(RedDist(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\RedDist.png')
        RedWS = watershed(RedDist);%splitting and indexing %it(RedWS(:,:,3))
        %imwrite(imadjust(uint16(max(RedWS,[],3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\RedWS.png')
        %%imwrite(imadjust(uint16(RedWS(:,:,3) > 0)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\RedWS.png')
        %imwrite(imadjust(logical(RedWS(:,:,3) > 0)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\RedWS.png')
        RedStencil = uint16(RedWS) .* uint16(RedBigMask); % Keep Big mask indexed pixels vol(RedStencil)
        RedFinalMask = RedStencil > 0; %vi(uint16(flip(RedFinalMask,1)))
        RedFinalMask = RedFinalMask | RedallMask; % Bring pixels from RedallMask back into RedFinalMask
        RedFinalMask = bwareaopen(RedFinalMask, SizeExclusionThreshold);
        RedallMask = RedFinalMask; %vi(uint16(flip(RedFinalMask,1))) vol(RedFinalMask)
        %imwrite(imadjust(uint16(max(RedallMask,[],3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\FinalDSredMask.png')
        %%imwrite(imadjust(uint16(RedallMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\FinalDSredMask.png')
        %% Segment high pH objects, mainly phagophores


        GreenAllDog1 = imfilter(double(ImG), fspecial('gaussian', 100, 1), 'symmetric') - imfilter(double(ImG), fspecial('gaussian', 100, 5), 'symmetric');%vi(flip(GreenAllDog1,1))
        %%imwrite(imadjust(uint16(GreenAllDog1(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\GreenDoG.png')
        GreenAllEdge = imfilter(double(ImG), fspecial('log', 20, 1), 'symmetric');%vi(flip(GreenAllEdge,1))
        %%imwrite(imadjust(im2uint16(GreenAllEdge(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\GreenEdge.png')

        GreenAllDog1Mask = GreenAllDog1 > 1000;%vi(double(flip(GreenAllDog1Mask,1)))
        %%imwrite(imadjust(uint16(GreenAllDog1Mask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\GreenDoGMask.png')
        GreenAllEdgeMask = bwareaopen((GreenAllEdge < -2000), 10);%vi(double(flip(GreenAllEdgeMask,1)))
        %imwrite(imadjust(uint16(GreenAllEdgeMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\GreenEdgeMask.png')
        %%imwrite(imadjust(uint16(GreenAllEdgeMask(:,:,3)) .* 2^16), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\GreenEdgeMask.png')

        GreenallMask3D = GreenAllDog1Mask | GreenAllEdgeMask; %vi(flip(uint8(GreenallMask3D),1))
        GreenallMask3D = bwareaopen(GreenallMask3D,10); %vi(flip(uint8(GreenallMask3D),1))
        %%imwrite(imadjust(uint16(GreenallMask3D(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\pHluorinMask.png')


        %GreenallMask = max(GreenallMask, [], 3);%vi(flip(uint8(GreenallMask),1))

        %% Compute ratio to be used for classification

        RatioIm = double(imfilter(ch1, fspecial('gaussian', 5, 2), 'symmetric')) ./ double(imfilter(ch2, fspecial('gaussian', 5, 2), 'symmetric'));
        % substitute NaN and Inf by numbers
        RatioIm(isinf(RatioIm)) = quantile(RatioIm(:), 0.99);
        RatioIm(isnan(RatioIm)) = quantile(RatioIm(:), 0.01); %vi(double(flip(RatioIm,1)))
        %%imwrite(uint16(2^16 .* (RatioIm(:,:,3)) ./ max(max(RatioIm(:,:,3)))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\RatioIm.png')
        %imtool(uint16(2^16 .* (RatioIm(:,:,3)) ./ max(max(RatioIm(:,:,3)))));
        %toto = uint16(2^16 .* (RatioIm(:,:,3)) ./ max(max(RatioIm(:,:,3))))
        
        %% Low pH candidates
        RatioImComp = imcomplement(RatioIm);
       
        RatioImComp = imtophat(RatioImComp, strel('disk', 15));
        RedVesicles15Mask = RatioImComp > 1.5;
        %vol(RatioImComp)
        %vol(RedVesicles15Mask)
        
         AutoLysoMask = zeros(size(ImR));
         
         RedVesicles15Objects = regionprops('table', RedVesicles15Mask, 'PixelIdxList');
         for v = 1:height(RedVesicles15Objects)
                 if iscell(RedVesicles15Objects{v,'PixelIdxList'}) 
                     CenterVec = ImG(RedVesicles15Objects{v,'PixelIdxList'}{:});
                 else
                     CenterVec = ImG(RedVesicles15Objects{v,'PixelIdxList'});
                 end
                 CenterMedian = median(CenterVec);
                 CenterMask = f_Create_Mask_from_ObjectList_Pixel_IDX(RedVesicles15Objects(v,:), 'PixelIdxList', ImG);%vol(CenterMask)
                 ContourMask = imdilate(CenterMask, strel('disk', 7)) - CenterMask; %vol(ContourMask)
                 if sum(ContourMask(:)) == 0
                     continue
                 end
                 ContourVec = ImG(logical(ContourMask));
                 ContourQuant = quantile(ContourVec, 0.5);
                 if ContourQuant > 1.5 * CenterMedian
                     AutoLysoMask(find(CenterMask)) = 1;
                 end
         end
         
         AutoLysoMask = bwareaopen(AutoLysoMask, 100); % vol(AutoLysoMask)
         
         
         
         %vol( AutoLysoMask)
         
         %vol(ImG + uint16(2^16*bwperim(RedVesicles15Mask)) + uint16(2^16*bwperim(RedVesicles15MaskDil)), 0, 5000)
         %vol(ImG + uint16(2^16*bwperim(AutoLysoMask)), 0, 5000)
        
        %% Fourrier and Euler
        for p = 1:5
            EulerImDoG(:,:,p) = f_LPF_by_FFT(GreenAllDog1(:,:,p), 'Butterworth', [10, 5], 0); % fourrier on DoG disconnects vesicle candidates from other contours
        end
        %%imwrite(imadjust(uint16(EulerImDoG(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\GreenDioG_FTB.png')
        %vi(double(flip(EulerImDoG,1)))

        AutophagosomeFourrierMask = EulerFiltering(EulerImDoG); %vi(double(flip(AutophagosomeFourrierMask,1)))
        %%imwrite(imadjust(uint16(AutophagosomeFourrierMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\FourrierEulerMask.png')

        %% Block of Hough transform approach for autophogsome filtering

        % Highlight rings for the hough transform
        GreenAllDog2 = imfilter(double(ch1), fspecial('gaussian', 100, 1), 'symmetric') - imfilter(double(ch1), fspecial('gaussian', 100, 25), 'symmetric');%vi(flip(GreenAllDog2,1))
        % Squeeze data into range [0 1]
        GreenAllDogNorm = GreenAllDog2 - min(GreenAllDog2(:)); % vi(flip(GreenAllDogNorm,1))
        GreenAllDogNormDouble = GreenAllDogNorm/max(GreenAllDogNorm(:));
        % vi(double(flip(GreenAllDogNormDouble,1)))

        GreenAllDogNormDoubleLP = imfilter(GreenAllDogNormDouble, fspecial('gaussian', 200, 20), 'symmetric'); %vi(double(flip(GreenAllDogNormDoubleLP,1)))
        GreenAllDogNormDoubleCorr = GreenAllDogNormDouble;
        % substitute bright phagophore pixels by local 'LP' gaussian weighted image average
        GreenAllDogNormDoubleCorr(GreenallMask3D) = GreenAllDogNormDoubleLP((GreenallMask3D)); 
        GreenAllDogNormDoubleCorr = imerode(GreenAllDogNormDoubleCorr, strel('disk',2));%vi(double(flip(GreenAllDogNormDoubleCorr,1)))

        % Apply hough transform to detect circles on each plane
        for p = 1:size(GreenAllDogNormDoubleCorr, 3)
            Autophagosome_Im = GreenAllDogNormDoubleCorr(:,:,p);
            [Autophagosome_centers, Autophagosome_radii, Autophagosome_metric] = imfindcircles(Autophagosome_Im, [3, 30], 'EdgeThreshold', 0.1, 'Sensitivity', 1);%Sensitivity:higher value->more circles detected;Edge:lower value->more circles found
            AutophagosomeCandidates = find(Autophagosome_metric > 0.2); % Keep candidates with good score only
            AutophagosomeMask(:,:,p) = f_Mask_DiscsfromHough(Autophagosome_Im, Autophagosome_centers, Autophagosome_radii, AutophagosomeCandidates); %vi(flip(uint8(A_Mask),1));
        end

        %imwrite(imadjust(uint16(AutophagosomeMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\HoughMask.png')
        AutophagosomeMask2D = sum(AutophagosomeMask,3); % check if object was found in multiple planes
        AutophagosomeMask2D = AutophagosomeMask2D >= 1;%vi(flip(uint8(AutophagosomeMask2D),1));
        AutophagosomeMask2D = cat(3, AutophagosomeMask2D,AutophagosomeMask2D,AutophagosomeMask2D,AutophagosomeMask2D,AutophagosomeMask2D);

        % Keep only those planes of the connected components which were detected via the Hough transform
        AutophagosomeMask = AutophagosomeMask .* AutophagosomeMask2D;
        %vi(flip(uint8(AutophagosomeMask),1));


        %% Autophagosome filtering for validation from Hough transform only (Fourrier/Euler are still excluded in this block)

        % Autophagosomes are round have no red bright spots but high median
        % absolute pixel deviation in the green channel i.e., the ring

        Autophagosome_CC = bwconncomp(AutophagosomeMask);
        Autophagosome_CandidatesConfirmed = regionprops('table', Autophagosome_CC, ch2, {'PixelIdxList','PixelValues'});
        Autophagosome_CandidatesConfirmedG = regionprops('table', Autophagosome_CC, ch1, {'PixelValues'});
        Autophagosome_CandidatesConfirmedG.Properties.VariableNames = {'pHluorinValues'};
        GreenMad = rowfun(@(x) mad(double(cell2mat(x)),1), Autophagosome_CandidatesConfirmedG, 'InputVariables', {'pHluorinValues'});% median absolute deviation
        GreenMad.Properties.VariableNames = {'GreenMad'};

        RedHigh95 = rowfun(@(x) quantile(x{:},0.9), Autophagosome_CandidatesConfirmed, 'InputVariables', {'PixelValues'});
        RedHigh95.Properties.VariableNames = {'RedHigh95'};
        Autophagosome_CandidatesConfirmed = [Autophagosome_CandidatesConfirmed, GreenMad, RedHigh95];

        % Autophagosome_CandidatesConfirmed = Autophagosome_CandidatesConfirmed(Autophagosome_CandidatesConfirmed.RedHigh95 < 110,:);
        Autophagosome_CandidatesConfirmed = Autophagosome_CandidatesConfirmed(Autophagosome_CandidatesConfirmed.RedHigh95 < 300,:);
        Autophagosome_CandidatesConfirmed = Autophagosome_CandidatesConfirmed(Autophagosome_CandidatesConfirmed.GreenMad > 20,:);
        AutophagosomeMask = f_Create_Mask_from_ObjectList_Pixel_IDX(Autophagosome_CandidatesConfirmed, 'PixelIdxList', AutophagosomeMask);%vi(flip(double(A_Mask2D) .* A_Mask .* double(RatioIm),1));
        AutophagosomeHoughMask = AutophagosomeMask;%vi(flip(uint8(AutophagosomeHoughMask),1));
        %vi(flip(uint8(AutophagosomeHoughMask + AutophagosomeMask),1));
        %% Autophagosome filtering for validation from pooled Hough transform and Fourrier/Euler

        AutophagosomeMask = AutophagosomeHoughMask | AutophagosomeFourrierMask;

            % Size filtering
            AutophagosomeMask = bwareaopen(AutophagosomeMask, 50);
            AutophagosomeMask = f_RemoveBigObjects(AutophagosomeMask, 10000); % vi(flip(uint8(AutophagosomeMask),1)) % do not go below 5000 

                %Shape filtering
                AutophagosomeMask = f_RemoveNonSphericalObjects(AutophagosomeMask, 1, 1.5); % vi(flip(double(AutophagosomeMask) .* double(ImR),1))
                AutophagosomeMaskDil = bwmorph(max(AutophagosomeMask,[],3), 'thicken', 9); % Dilate but keep proximate autophagosomes split
                AutophagosomeMaskDil = cat(3, AutophagosomeMaskDil, AutophagosomeMaskDil, AutophagosomeMaskDil, AutophagosomeMaskDil, AutophagosomeMaskDil); % vi(uint8(flip(AutophagosomeMaskDil,1)))

                    % Evaluate the ratio intensity between vesicle border and vesicle center
                    AutophagosomeMask_Centers = imerode(AutophagosomeMaskDil, strel('disk', 5));
                    AutophagosomeMask_Rings = AutophagosomeMaskDil - AutophagosomeMask_Centers;

                    AutophagosomeCC = bwconncomp(AutophagosomeMaskDil,6);
                    AutophagosomeStencil = bwlabeln(AutophagosomeMaskDil,6);%vi(flip(AutophagosomeStencil,1));
                    ObjectsAutophagosomeLabels = regionprops('table', AutophagosomeCC, AutophagosomeStencil, {'PixelIdxList','MaxIntensity'});

                    thisRatio = [];
                    decisionVec = [];

                    for v = 1:height(ObjectsAutophagosomeLabels)

                        thisObjectRow = ObjectsAutophagosomeLabels(v,:);
                        thisVesicle = f_Create_Mask_from_ObjectList_Pixel_IDX(thisObjectRow, 'PixelIdxList', AutophagosomeMaskDil); %vi(flip(thisVesicle,1))
                        thisRing = thisVesicle & AutophagosomeMask_Rings; %vi(flip(uint8(thisRing),1))
                        thisCenter = thisVesicle & AutophagosomeMask_Centers; %vi(flip(uint8(thisCenter),1))
                        thisRatio(v) = quantile(RatioIm(thisRing),0.75) / quantile(RatioIm(thisCenter), 0.75);% heuristic threshold is 1.3 20161103 Javier and Paul
                        if thisRatio(v) > 0.95
                            decisionVec(v) = 1;
                        else
                            decisionVec(v) = 0;
                        end

                    end

                    ObjectsAutophagosomeLabels = [ObjectsAutophagosomeLabels, array2table(decisionVec')];
                    ObjectsAutophagosomeLabels.Properties.VariableNames(end) = {'Decision'};
                    ObjectsAutophagosomeLabels = ObjectsAutophagosomeLabels(ObjectsAutophagosomeLabels.Decision == 1, :);
                    AutophagosomeMask = f_Create_Mask_from_ObjectList_Pixel_IDX(ObjectsAutophagosomeLabels, 'PixelIdxList', RatioIm); %vi(flip( RatioIm .* AutophagosomeMask, 1))


        % Restore 3D smoothness            
        AutophagosomeMask = AutophagosomeMask .* (AutophagosomeHoughMask | AutophagosomeFourrierMask);
        %vi(uint8(flip(AutophagosomeMask, 1)))
        %imwrite(imadjust(uint16(AutophagosomeMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\AutophagosomeMask.png')
        %%imwrite(imadjust(uint16(AutophagosomeMask(:,:,3)) .* 2^16), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\AutophagosomeMask.png')

        %%%%%%%%%%%%%%%%%%%%
        %% Classification %%
        %%%%%%%%%%%%%%%%%%%%

        % Prevent fusions between the AutophagosomeMask and other masks
        AutophagosomesPerim = f_bwperim3D(AutophagosomeMask); %vi(uint8(flip(AutophagosomesPerim,1))); vol(AutophagosomeMask)
        %vol(ImR + 2^16*uint16(AutophagosomesPerim), 0, 5000)
        
        RedThreshold = quantile(ImR(RedallMask), 0.5);

        AutophagosomeValidationMask = (AutophagosomeMask .* (ImR > RedThreshold));
        APvalObjects = regionprops('table', bwconncomp(AutophagosomeMask), AutophagosomeValidationMask, {'Area', 'PixelValues', 'PixelIdxList'});
        BrightPixels = rowfun(@(x) sum(x{:}), APvalObjects, 'InputVariables', 'PixelValues');
        BrightPixels.Properties.VariableNames = {'BrightPixels'};
        APvalObjects = [APvalObjects, BrightPixels];
        BrightRatio = rowfun(@(x,y) x/y, APvalObjects, 'InputVariables', {'BrightPixels', 'Area'});
        BrightRatio.Properties.VariableNames = {'BrightRatio'};
        APvalObjects = [APvalObjects, BrightRatio];
        APvalObjectsConfirmed = APvalObjects(APvalObjects.BrightRatio > 0.5,:);
        AutophagosomeMask = f_Create_Mask_from_ObjectList_Pixel_IDX( APvalObjectsConfirmed, 'PixelIdxList', ImR);
        %vol(AutophagosomeMask)
         
        RedallMaskPerim = f_bwperim3D(RedallMask);
        AutophagosomesProtector = ~imdilate(AutophagosomesPerim, strel('sphere',1));
        RedVesicleProtector = ~imdilate(RedallMaskPerim, strel('sphere',1));
        AutolysosomeProtector = ~imdilate(AutoLysoMask, strel('sphere',1)); %vol(AutolysosomeProtector)

        % Create a vesicle mask combining all red fluorescent vesicles and the autophagosomes
        %VesiclesAll = (AutophagosomeMask | RedallMask) & AutophagosomesProtector; %vi(uint8(flip(VesiclesAll,1)));
        VesiclesAll = ((AutophagosomeMask | RedallMask | AutoLysoMask) & AutophagosomesProtector) | ((AutophagosomeMask | RedallMask | AutoLysoMask) & AutolysosomeProtector); %vol(VesiclesAll);
        VesiclesAll = f_imclearborder3D(VesiclesAll);
        %%imwrite(imadjust(uint16(VesiclesAll(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\VesiclesAll.png')

        % Create a feature table for the classification
        VesiclesAllCC = bwconncomp(VesiclesAll); % toto = labelmatrix(VesiclesAllCC);VesicleObjects=VesicleObjects(210,:)%vol(toto)

        VesicleObjects0 = regionprops('table', VesiclesAllCC, ch2, {'PixelIdxList', 'MeanIntensity', 'Image'});
        VesicleObjects0 = VesicleObjects0(VesicleObjects0.MeanIntensity < 300,:);
        LowRedVesicleMask = f_Create_Mask_from_ObjectList_Pixel_IDX(VesicleObjects0, 'PixelIdxList', ch2);
        VesiclesAll = VesiclesAll & ~LowRedVesicleMask; % vol(VesiclesAll)
        VesiclesAllCC = bwconncomp(VesiclesAll); % Overwrite
        VesicleObjects0 = regionprops('table', VesiclesAllCC, ch2, {'PixelIdxList', 'MeanIntensity', 'Image'});
        
        ObjectsPlusEccentricity = VesicleEccentricity(VesicleObjects0);
        ObjectsPlusEccentricity = ObjectsPlusEccentricity(ObjectsPlusEccentricity.Eccentricity>0.9,:); % Non circular vesicles
        NonCircularVesicleMask = f_Create_Mask_from_ObjectList_Pixel_IDX(ObjectsPlusEccentricity, 'PixelIdxList', ch2);
        %imwrite(imadjust(im2uint16(NonCircularVesicleMask(:,:,3))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\NonCircularVesicleMask.png')
        %vol(MaskBad)
        %vol(NonCircularVesicleMask)
        
        %vol(ch2)
        VesicleObjects1 = regionprops('table', VesiclesAllCC, RedallMask, {'PixelIdxList', 'PixelValues', 'Image'});

        
        VesicleObjects1.Properties.VariableNames{end} = 'redMaskValues';
        VesicleObjects2 = regionprops('table', VesiclesAllCC, GreenallMask3D, {'PixelValues'});
        VesicleObjects2.Properties.VariableNames{end} = 'greenMaskValues';
        VesicleObjects3 = regionprops('table', VesiclesAllCC, AutophagosomeMask, {'PixelValues'});
        VesicleObjects3.Properties.VariableNames{end} = 'autophagosomeMaskValues';
        VesicleObjects4 = regionprops('table', VesiclesAllCC, NonCircularVesicleMask, {'PixelValues'});
        VesicleObjects4.Properties.VariableNames{end} = 'NonCircularMaskValues';
        VesicleObjects5 = regionprops('table', VesiclesAllCC, RatioIm, {'PixelValues'});
        VesicleObjects5.Properties.VariableNames{end} = 'RatioValue';
        VesicleObjects6 = regionprops('table', VesiclesAllCC, ImG, {'PixelValues'});
        VesicleObjects6.Properties.VariableNames{end} = 'ImG';   
        VesicleObjects7 = regionprops('table', VesiclesAllCC, ImR, {'PixelValues'});
        VesicleObjects7.Properties.VariableNames{end} = 'ImR';   
        VesicleObjects8 = regionprops('table', VesiclesAllCC, AutoLysoMask, {'PixelValues'});
        VesicleObjects8.Properties.VariableNames{end} = 'AutoLysoMask';
        
        VesicleObjects9 = regionprops('table', VesiclesAllCC, {'SubArrayIdx'});
        VesicleObjects9.Properties.VariableNames(:) = {'SubIndexes'};
        
        VesicleObjects10 = regionprops('table', VesiclesAllCC, {'Image'});
        
        
        
        %% Get Padded Subindexes
        XRange = rowfun(@(x) {x{1}}, VesicleObjects9);
        YRange = rowfun(@(x) {x{2}}, VesicleObjects9);
        ZRange = rowfun(@(x) {x{3}}, VesicleObjects9);
        
        PaddingRange = 1;
        PerimeterRange = 1;
        XRange = rowfun(@(x) {[(min(x{:})-PaddingRange) : (max(x{:})+PaddingRange)]}, XRange); XRange.Properties.VariableNames{1} = 'PosX';
        YRange = rowfun(@(x) {[(min(x{:})-PaddingRange) :  (max(x{:})+PaddingRange)]}, YRange); YRange.Properties.VariableNames{1} = 'PosY';        
        %ZRange = rowfun(@(x) {[min(x{:})-0, x{:}, max(x{:})+0]}, ZRange); ZRange.Properties.VariableNames{1} = 'Z';
        
        SubIdxPadded = [XRange, YRange, ZRange];
        PaddedImages = rowfun(@(x) {padarray(x{:}, [PaddingRange, PaddingRange, 0])}, VesicleObjects10, 'InputVariables', {'Image'});
        PaddedImages.Properties.VariableNames{1} = 'ImagePadded';
        PaddedImagesDilated = rowfun(@(x) {imdilate(x{:}, strel('disk',PaddingRange))}, PaddedImages, 'InputVariables', {'ImagePadded'});
        PaddedImagesDilated.Properties.VariableNames{1} = 'ImagePaddedDilated';
        PaddedImagesEroded =  rowfun(@(x) {imerode(x{:}, strel('disk',PerimeterRange))}, PaddedImages, 'InputVariables', {'ImagePadded'});
        PaddedImagesEroded.Properties.VariableNames{1} = 'ImagePaddedEroded';
        %vol(PaddedImages{3,1}{1,1})
        %vol(PaddedImagesDilated{3,1}{1,1})
        %vol(PaddedImagesEroded{3,1}{1,1})
        SingleVesicleMasks = [PaddedImages, PaddedImagesDilated, PaddedImagesEroded];
        PerimeterMasks = rowfun(@(x,y) {x{:}-y{:}}, SingleVesicleMasks, 'InputVariables', {'ImagePadded', 'ImagePaddedEroded'}); 
        PerimeterMasks.Properties.VariableNames{1} = 'Perimeter';
        %vol(PerimeterMasks{3,1}{1,1})
        
        OuterMasks = rowfun(@(x,y) {x{:}-y{:}}, SingleVesicleMasks, 'InputVariables', {'ImagePaddedDilated', 'ImagePadded'}); 
        OuterMasks.Properties.VariableNames{1} = 'Outside';
        %vol(OuterMasks{3,1}{1,1})
        
        %SingleVesicleGreenImBlocks = rowfun(@(x) {ImG(x{:})}, cell2table({SubIdxPadded{:,1}, SubIdxPadded{:,2}, SubIdxPadded{:,3}}));
        SingleVesicleGreenImBlocks = rowfun(@(a,b,c) {ImG([a{:}], [b{:}], [c{:}])}, SubIdxPadded, 'InputVariables', SubIdxPadded.Properties.VariableNames(:)');
        SingleVesicleGreenImBlocks.Properties.VariableNames{end} = 'Im';
        %vol(SingleVesicleGreenImBlocks{3,1}{1,1})
        
        SingleVesicleZones = [PaddedImagesEroded, PerimeterMasks, OuterMasks, SingleVesicleGreenImBlocks];
        GreenPerim = rowfun(@(x,y) getMedianWithinMask(x,y), SingleVesicleZones, 'InputVariables', {'Im','Perimeter'});
        GreenPerim.Properties.VariableNames{end} = 'PerimMedian';
        GreenOutside = rowfun(@(x,y) getMedianWithinMask(x,y), SingleVesicleZones, 'InputVariables', {'Im','Outside'});
        GreenOutside.Properties.VariableNames{end} = 'OutsideMedian';
        GreenCenter = rowfun(@(x,y) getMedianWithinMask(x,y), SingleVesicleZones, 'InputVariables', {'Im','ImagePaddedEroded'});
        GreenCenter.Properties.VariableNames{end} = 'CenterMedian';

        %VesicleObjects9{1,:}
        
        %vol(toto{1,1}{1,1})
        
        
        
        %VesicleObjects = [VesicleObjects1, VesicleObjects2, VesicleObjects3, VesicleObjects4, VesicleObjects5, VesicleObjects6, VesicleObjects7, VesicleObjects8];
        VesicleObjects = [VesicleObjects1, VesicleObjects2, VesicleObjects3, VesicleObjects4, VesicleObjects5, VesicleObjects6, VesicleObjects7, VesicleObjects8, GreenPerim, GreenOutside, GreenCenter];
  
        %%%%%%%%%%%%%%%%%%%%%%
        %Classify phagophores or anything else %%%%% 2 classes %%%%%%
        
        LowpHDecision = rowfun(@(x,y,z) LowpHDecider(x,y,z), VesicleObjects, 'InputVariables', {'autophagosomeMaskValues', 'NonCircularMaskValues', 'RatioValue'});
        LowpHDecision.Properties.VariableNames = {'LowpHDecision'};
        VesicleObjects = [VesicleObjects, LowpHDecision];
        LowpHVesicleMask = VesicleObjects(VesicleObjects.LowpHDecision == 1,:);
        LowpHVesicleMask = f_Create_Mask_from_ObjectList_Pixel_IDX(LowpHVesicleMask, 'PixelIdxList', ch2);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%
        % Take classification decisions
        AutophagosomeDecision = rowfun(@(x,y) AutoPhagosomeDecider(x,y), VesicleObjects, 'InputVariables', {'autophagosomeMaskValues', 'NonCircularMaskValues'});
        AutophagosomeDecision.Properties.VariableNames = {'AutophagosomeDecision'};
        VesicleObjects = [VesicleObjects, AutophagosomeDecision];

        PhagophoreDecision = rowfun(@(x,y,z,omega,ImG,ImR,AutoLysoMask,CenterGreen,OutsideGreen) PhagophoreDecider(x,y,z,omega,ImG,ImR,AutoLysoMask,CenterGreen, OutsideGreen), VesicleObjects, 'InputVariables', {'redMaskValues', 'greenMaskValues', 'AutophagosomeDecision', 'RatioValue', 'ImG', 'ImR', 'AutoLysoMask', 'CenterMedian', 'OutsideMedian'});
        PhagophoreDecision.Properties.VariableNames = {'PhagophoreDecision'};
        VesicleObjects = [VesicleObjects, PhagophoreDecision];

        LateAutolysosomeDecision = rowfun(@(x,y,z,l,omega,CenterGreen,PerimGreen) LateAutoysosomeDecider(x,y,z,l,omega,CenterGreen,PerimGreen), VesicleObjects, 'InputVariables', {'redMaskValues', 'greenMaskValues', 'AutophagosomeDecision', 'PhagophoreDecision', 'RatioValue', 'CenterMedian', 'PerimMedian'});
        LateAutolysosomeDecision.Properties.VariableNames = {'LateAutolysosomeDecision'};
        VesicleObjects = [VesicleObjects, LateAutolysosomeDecision];

        EarlyAutolysosomeDecision = rowfun(@(x,y,z,l,m,AutoLysoMask,CenterGreen,PerimGreen) EarlyAutolysosomeDecider(x,y,z,l,m,AutoLysoMask,CenterGreen,PerimGreen), VesicleObjects, 'InputVariables', {'redMaskValues', 'greenMaskValues', 'AutophagosomeDecision', 'PhagophoreDecision', 'LateAutolysosomeDecision', 'AutoLysoMask', 'CenterMedian', 'PerimMedian'});
        EarlyAutolysosomeDecision.Properties.VariableNames = {'EarlyAutolysosomeDecision'};
        VesicleObjects = [VesicleObjects, EarlyAutolysosomeDecision];

        % Test if unclassified vesicles are remaining
        GroupSum = rowfun(@(a,b,c,d) a+b+c+d, VesicleObjects, 'InputVariables', {'AutophagosomeDecision', 'PhagophoreDecision', 'LateAutolysosomeDecision', 'EarlyAutolysosomeDecision'});
        GroupSum.Properties.VariableNames = {'GroupSum'};
        VesicleObjects = [VesicleObjects, GroupSum];
     
        % Analyse vesicle size
        Diameters = rowfun(@(r) GetDiameter(r), VesicleObjects1,'InputVariables', {'Image'});
        VesicleObjects = [VesicleObjects, Diameters];
        VesicleObjects.Properties.VariableNames{end} = 'MajorAxisLength';
        
        BarCodeThis = unique(InfoTable.Barcode); BarCodeThis = BarCodeThis{:};
        ColumnThis = InfoTable.Column(i);
        RowThis = InfoTable.Row(i);
        FieldThis = InfoTable.field(i);
        AreaNameThis = InfoTable.AreaName{i};
        
        VesicleObjectsOutput = VesicleObjects(:, 7:end);
        CellArrayDummy = cell(height(VesicleObjectsOutput), 1);
        BarCodeColumn = CellArrayDummy;
        BarCodeColumn = cellfun(@(x) BarCodeThis, BarCodeColumn, 'UniformOutput', false);
        AreaNameColumn = CellArrayDummy;
        AreaNameColumn = cellfun(@(x) AreaNameThis, AreaNameColumn, 'UniformOutput', false);        
        ColumnColumn = CellArrayDummy;
        ColumnColumn = cellfun(@(x) ColumnThis, ColumnColumn, 'UniformOutput', false);
        RowColumn = CellArrayDummy;
        RowColumn = cellfun(@(x) RowThis, RowColumn, 'UniformOutput', false);
        FieldColumn = CellArrayDummy;
        FieldColumn = cellfun(@(x) FieldThis, FieldColumn, 'UniformOutput', false);
        VesicleObjectsOutput = [VesicleObjectsOutput, BarCodeColumn, AreaNameColumn, ColumnColumn, RowColumn, FieldColumn];
        VesicleObjectsOutput.Properties.VariableNames(15:end) = {'Barcode','AreaName','Column','Row','Field'};
        
        %% Append tables while looping of barcodes and their images
        if progress == 1
            ObjectsAll = VesicleObjectsOutput;
        else
            ObjectsAll = [ObjectsAll; VesicleObjectsOutput];
        end

        % Display previews showing classified vesicles
        if showIllustration
            [LabelMask] = LabelCategories(VesicleObjects, ImG); % vol(LabelMask)
            [LabelMask] = LabelCategories(VesicleObjects, ImR);
        end
        
        
        %% Previews
        
         % Scalebar
         imSize = [size(ch1,1),size(ch1,2)];
        [BarMask, BarCenter] = f_barMask(20, 0.2152, imSize, imSize(1)-50, 50, 5); %it(BarMask)
        
        [AutophagosomeMask, PhagophoreMask, LateAutolysosomeMask, EarlyAutolysosomeMask] = LabelMasks(VesicleObjects, ch1);
    
        Preview1 = imoverlay2(f_ImAdjust(max(ch1, [], 3)), bwperim(max(AutophagosomeMask,[],3)), [1 0 0]); % pHluorin
        Preview1 = imoverlay2(Preview1, bwperim(max(PhagophoreMask,[],3)), [0 1 0]); % pHluorin
        Preview1 = imoverlay2(Preview1, bwperim(max(LateAutolysosomeMask,[],3)), [0 0 1]); % pHluorin
        Preview1 = imoverlay2(Preview1, bwperim(max(EarlyAutolysosomeMask,[],3)), [1 1 0]); % pHluorin
        Preview1 = imoverlay2(Preview1, BarMask, [1 1 1]); % pHluorin
        %it(Preview1)
        %imwrite(Preview1, 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\ClassesOnPHluorin.png')
        

        Preview2 = imoverlay2(f_ImAdjust(max(ch2, [], 3), 0.985), bwperim(max(AutophagosomeMask,[],3)), [1 0 0]); % pHluorin
        Preview2 = imoverlay2(Preview2, bwperim(max(PhagophoreMask,[],3)), [0 1 0]); % pHluorin
        Preview2 = imoverlay2(Preview2, bwperim(max(LateAutolysosomeMask,[],3)), [0 0 1]); % pHluorin
        Preview2 = imoverlay2(Preview2, bwperim(max(EarlyAutolysosomeMask,[],3)), [1 1 0]); % pHluorin
        Preview2 = imoverlay2(Preview2, BarMask, [1 1 1]); % pHluorin
        %it(Preview2)


        
        PreviewRGB = cat(3, imadjust(ch2(:,:,3),[0 0.02]), imadjust(ch1(:,:,3),[0 0.06]), zeros([size(ch1,1), size(ch1,2)], 'uint16'));
        %it(PreviewRGB)
        %imwrite(PreviewRGB, 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\RGB.png')
        
        PreviewPath1 = [SavePath filesep 'Previews' filesep num2str(RowThis) '_' num2str(ColumnThis) '_' num2str(FieldThis) '_' AreaNameThis '_ch1'];
        PreviewPath2 = [SavePath filesep 'Previews' filesep num2str(RowThis) '_' num2str(ColumnThis) '_' num2str(FieldThis) '_' AreaNameThis '_ch2'];
        PreviewPathRGB = [SavePath filesep 'Previews' filesep num2str(RowThis) '_' num2str(ColumnThis) '_' num2str(FieldThis) '_' AreaNameThis '_RGB'];
        
        imwrite(Preview1, [PreviewPath1, '.png'])
        imwrite(Preview2, [PreviewPath2, '.png'])
        imwrite(PreviewRGB, [PreviewPathRGB, '.png'])
        
    end % images
end % Barcodes

writetable(ObjectsAll(:,4:end), [SavePath, filesep, 'ObjectsAll.csv'])
save([SavePath, filesep, 'ObjectsAll.mat'], 'ObjectsAll')


