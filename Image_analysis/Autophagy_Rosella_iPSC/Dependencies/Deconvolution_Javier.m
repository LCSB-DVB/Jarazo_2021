function [ImR, ImG] = Deconvolution_Javier(CreateNewPSF, ch1, ch3)
%We take each channel separate to deconvolve them with a Richards and Wolf
%psf Model. The modelling of the psf was done once and the parameters saved
%as a screenshot and as an psf.mat This model is accounting for several
%properties of the microscope but not the confocallity. 
 


%% Deconvolution

if CreateNewPSF == 1
    
    % Manual setup of psfs at project start
    % Go to http://bigwww.epfl.ch/algorithms/psfgenerator/
    % download PSFGenerator.jar and add it to your workspace
    javaaddpath ./PSFGenerator.jar % add java library for the generation of model based point spread functions to matlab's java path
    PSFGenerator.gui; % Set the parameters of the GUI. You may want to set the psf z size to the same size as your image Im's z size
    % See screenshot GUI_PSF_20161913.PNG
    psf = PSFGenerator.get; % Compute a PSF from the GUI parameters
    if decision3D == 1
        disp('3D psf in use')
    else
        psf = psf(:,:,2);
    end
    psf = double(psf);
    %vi(psf)
    save('psf.mat', 'psf')% Please rename manually to indicate channel and date. PS documentatation by snipping tool is recommended
    
else
    
    psf_DR = load('psfDeepRed_5planesBin2.mat');
    psf_DR = psf_DR.psf;
    
    psf_Red = load('psfRed_5planesBin2.mat');
    psf_Red = psf_Red.psf;
    
    psf_Green = load('psfGreen_5planesBin2.mat');
    psf_Green = psf_Green.psf;
    
end

ImR = double(ch1); % image data type should match psf data type
[RedDeconvolvedIm,PSF] = deconvblind(ImR, psf_Red, 10, 0); % Matlab internal function
%RedDeconvolvedIm = uint16(2^16 .* RedDeconvolvedIm ./ (max(RedDeconvolvedIm(:))));
ImR = uint16(2^16 .* RedDeconvolvedIm ./ (max(RedDeconvolvedIm(:))));

ImG = double(ch3); % image data type should match psf data type
[GreenDeconvolvedIm,PSF] = deconvblind(ImG, psf_Green, 10, 0); % Matlab internal function
%GreenDeconvolvedIm = uint16(2^16 .* GreenDeconvolvedIm ./ (max(GreenDeconvolvedIm(:))));
ImG = uint16(2^16 .* GreenDeconvolvedIm ./ (max(GreenDeconvolvedIm(:))));
%vi(ImG)

end

