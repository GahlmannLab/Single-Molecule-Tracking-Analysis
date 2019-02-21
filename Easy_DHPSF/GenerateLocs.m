function [PSFLocs,numPSFLocs,maxPeakImg] = GenerateLocs(cropHeight, cropWidth, gaussianFilter,dataFT,templateFT,FOVmask,peakThreshold,numTemplates,minDistBetweenSMs,nhaData,trueFOVmask, FTFOVmask)

PSFLocs = [];
maxPeakImg = zeros(cropHeight,cropWidth);
% matrix PSFLocs stores information about double helices that were
% found via template matching
% rows are different matches
% [xLocation yLocation matchingTemplateNumber matchConfidence];
tempPSFLoc = cell(numTemplates,1);
numPSFLocs = 0;

for b=1:numTemplates
    % try no prefiltering
    %H = 1;
    % try phase correlation
    %H = 1./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
    % try weighted phase correlation (emphasizing low frequency
    % components
    if nhaData
        H = gaussianFilter;
    else
        %if dataFT has zeros, make small number to avoid dividing by 0
        if sum(sum(dataFT == 0)) ~=0
            [idx1 idx2] = find(dataFT == 0);
            dataFT(idx1, idx2) = 1/1000000;
        end
        H = gaussianFilter./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
    end
    % normalize H so it doesn't add any energy to template match
    %H = H / sqrt(sum(abs(H(:)).^2));
    
    peakImg = ifftshift(ifft2(dataFT.*squeeze(templateFT(b,:,:)).*H));
    if trueFOVmask == 1 % prevents finding matches outside FOV (due to e.g. wrapping around)
        %if exist('FOVmask')
        peakImg=peakImg.*FOVmask;
    end
    % normalize response of peakImg by dividing by number of pixels in
    % data
    %peakImg = peakImg / (cropHeight*cropWidth);
    maxPeakImg = max(maxPeakImg, peakImg);
    
    %threshold = mean(peakImg(:))+peakThreshold*std(peakImg(:));
    peakImg(peakImg < peakThreshold(:,b)) = peakThreshold(:,b);
    peakImg(~FTFOVmask) = peakThreshold(:,b);%This mask removes the excess brightness introduced by the fft/ifft of the FOV
    if isreal(peakImg) && sum(sum(isnan(peakImg)))==0
        temp = find(imregionalmax(peakImg));
    else
        peakImg(isnan(peakImg)) = 0; %inserted to deal with NaNs -AC 6/22
        temp = find(imregionalmax(real(peakImg)));
    end
    
    
    % make sure threshold didn't eliminate all peaks and create
    % lots of matches
    if length(temp) < cropHeight*cropWidth/2;
        [tempY, tempX] = ind2sub([cropHeight cropWidth],temp);
        tempPSFLoc{b,1} = [tempX tempY b*ones(length(temp),1) peakImg(temp)];
        numPSFLocs = numPSFLocs+length(temp);
    end
end
%clear H dataFT peakImg
%% filter out extraneous matches due to very strong signals
if numPSFLocs > 0
    % sort location matrix in decending order of confidence
    tempPSFLocs = cell2mat(tempPSFLoc);
    temp = sortrows(tempPSFLocs,-4);
    % copy most confident match to list of locations
%     tempDist = bsxfun(@minus,temp(:,1),temp(:,1)').^2+bsxfun(@minus,temp(:,2),temp(:,2)').^2;
%     tempDist2 = tempDist < minDistBetweenSMs^2;
%     [~,temp1] = max(tempDist2);
% %     find(tempDist2,1)
%     PSFLocsTK = temp(temp1 == (1:length(temp1)),:);
    tempPSFLocs(1,:) = temp(1,:);
    numPSFLocs = 1;
    for b=2:size(temp,1)
        % make sure that this candidate location is a minimum distance away
        % from all other candidate locations
        %if sum((temp(b,1)-PSFLocs(1:numPSFLocs,1)).^2 + (temp(b,2)-PSFLocs(1:numPSFLocs,2)).^2 >= minDistBetweenSMs^2) == numPSFLocs
        if sum((temp(b,1)-tempPSFLocs(1:numPSFLocs,1)).^2 + (temp(b,2)-tempPSFLocs(1:numPSFLocs,2)).^2 >= minDistBetweenSMs^2) == numPSFLocs
            % add it to list of locations
            numPSFLocs = numPSFLocs + 1;
            tempPSFLocs(numPSFLocs,:) = temp(b,:);
        end
    end
PSFLocs = tempPSFLocs(1:numPSFLocs,:);

end
end