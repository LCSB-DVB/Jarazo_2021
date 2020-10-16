function [HPFim] = f_HPF_by_FFT(Im, FilterType, FilterParameters, illustrate)
%Perform fourrier transform and butterworth high pass filtering
%Author: Paul Antony 2016/10/31
%   Im: the input image
%   FilterType: currently only 'Butterworth'
%   FilterParameters: row vector with two elements:
%   1) cutoff frequeny
%   2) order of the filter
%   illustrate: show illustrations (1) or not (0)
%   For more details on fourrier transforms see http://www.ogemarques.com/books#ffs-tabbed-14
%   example of use: HPF = f_LPF_by_FFT(MyGreatImage, 'Butterworth', [50, 1], 1);

[M,N] = size(Im);
dist = distmatrix(M,N);

ImD = im2double(Im);
ImD_fft = fft2(ImD);

if strcmp(FilterType, 'Butterworth')
    cutoff = FilterParameters(1); % cutoff frequency
    order =  FilterParameters(2);% order of the filter
    H_but = 1 ./ (1 + (cutoff ./ dist) .^ (2*order));
    ImFilt = H_but .* ImD_fft;
    HPFim = real(ifft2(ImFilt));

    if illustrate == 1
        figure;
        subplot(2,3,1), imshow(Im,[]); title('Input image')
        subplot(2,3,2), imshow(fftshift(log(1 + abs(ImD_fft))),[]); title('Fourrier transform')
        subplot(2,3,3), imshow(fftshift(log(1 + H_but)), []); title([FilterType, ' filter'])
        subplot(2,3,4), imshow(fftshift(log(1 + abs(ImFilt))), []); title('Filtered Fourrier transform')
        subplot(2,3,5), imshow(HPFim, []); title('Low pass filtered image')
        subplot(2,3,6), imshow(Im+HPFim, [], 'Colormap', gray); title('Enhanced image')
        annotation('textbox',[0.1 0.9 0.8 0.1],'String','Low pass filtering by fourrier transform and butterworth filtering', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'LineStyle', 'none')
        if 0
            imwrite(imadjust(uint16(fftshift(log(1 + abs(ImD_fft))))), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\FourrierTransform.png')
            imwrite(imadjust(im2uint16(fftshift(log(1 + H_but))),[]), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\Butterworth.png')
            imwrite(imadjust(uint16(HPFim)), 'S:\HCS_Platform\PaulAntony\AutophagyStaging\LC3_Staging\FourrierIm.png')
        end
    end
       
end

end

