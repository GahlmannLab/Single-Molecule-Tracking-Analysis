function PlotImage(PSFfits,PSFLocs,cropWidth,cropHeight,ROI,data,bkgndImg,dispRaw,printOutputFrames,numPSFLocs,peakThreshold,maxPeakImg,templateColors,frame,stack,hSMFits,iROI)


% if the fit is good, add it to the reconstructed image
% okayFit = PSFfits(:,11) > 0;
[xIdx, yIdx] = meshgrid(1:cropWidth,1:cropHeight);
reconstructImg = zeros(cropHeight, cropWidth);
PSFfits(:,3) = PSFfits(:,3)- ROI(1)+1-iROI(1)+1;
PSFfits(:,4) = PSFfits(:,4)- ROI(2)+1-iROI(2)+1;
PSFfits(:,5) = PSFfits(:,5)- ROI(1)+1-iROI(1)+1;
PSFfits(:,6) = PSFfits(:,6)- ROI(2)+1-iROI(2)+1;
for b = 1:size(PSFfits,1)
    if PSFfits(b,11) > 0
        reconstructImg = reconstructImg + ...
            PSFfits(b,1).*exp( -((xIdx-PSFfits(b,3)).^2+(yIdx-PSFfits(b,4)).^2.) / (2.*PSFfits(b,7).^2)) ...
            +PSFfits(b,2).*exp( -((xIdx-PSFfits(b,5)).^2+(yIdx-PSFfits(b,6)).^2.) / (2.*PSFfits(b,8).^2));
    end
end
%%  plot results of template matching and fitting
if dispRaw
    dataToPlot= data+bkgndImg;
    reconstToPlot = reconstructImg + bkgndImg;
    dataTitle = ['Frame ' num2str(frame) ': raw data - darkAvg counts'];
else
    dataToPlot = data;
    reconstToPlot = reconstructImg;
    dataTitle = ['Frame ' num2str(frame) ': raw data - darkAvg & bkgnd'];
end
set(0,'CurrentFigure',hSMFits);
subplot('Position',[0.025 0.025 .85/3 .95]);
imagesc(maxPeakImg,[0 3*min(peakThreshold(stack,:))]);axis image;
title({'Peaks correspond to likely template matches' ...
    [num2str(numPSFLocs) ' matches found']});

subplot('Position',[0.075+.85/3 0.025 .85/3 .95]);
imagesc(dataToPlot);axis image;colormap hot;
hold on;
for b=1:numPSFLocs
    if PSFfits(b,11) > 0
    %plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
    %    'MarkerSize', 15*PSFLocs(b,4)/peakThreshold(b), ...
    %    'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
    plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
        'MarkerSize', 15*PSFLocs(b,4)/peakThreshold(stack,PSFLocs(b,3)), ...
        'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
    %Changed peakThreshold(b) to peakThreshold(PSFLocs(b,3)) because b
    %does not seem to logically correspond to the correct template, whereas the third
    %column of PSFLocs is defined as corresponding to a specific
    %template. Furthermore, whenever numPSFLocs > length(peakThreshold)
    %there is an error. -AC 6/22
    end
end
hold off;
title({dataTitle ...
    ['ROI [xmin ymin width height] = ' mat2str(ROI)]});

subplot('Position',[0.125+2*.85/3 0.025 .85/3 .95]);
%         imagesc(reconstructImg+bkgndMean,[min(data(:)) max(data(:))]);axis image;
% imagesc(reconstToPlot,[min(dataToPlot(:)) max(dataToPlot(:))]);axis image;
imagesc(reconstToPlot);axis image;
title({'Image reconstructed from fitted matches' ...
    [num2str(sum(PSFfits(:,11)>0)) ' successful fits']});
drawnow;
if printOutputFrames == 1
    set(gcf,'PaperPositionMode','auto');
    saveas(hSMFits, ['output images' filesep 'frame ' num2str(frame) '.tif']);
end

PSFfits(:,3) = PSFfits(:,3) + ROI(1)-1+iROI(1)-1;
PSFfits(:,4) = PSFfits(:,4) + ROI(2)-1+iROI(2)-1;
PSFfits(:,5) = PSFfits(:,5) + ROI(1)-1+iROI(1)-1;
PSFfits(:,6) = PSFfits(:,6) + ROI(2)-1+iROI(2)-1;
end
