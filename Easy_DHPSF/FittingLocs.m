function PSFfits = FittingLocs(PSFLocs,FOVmask1,ROI,boxRadius,fittingMethod,cropWidth,cropHeight,templateLocs,templateSize,data,datavar,bkgndImg,sigmaBounds,conversionFactor,nmPerPixel,sigSearchBounds,sigmaRatioLimit,lobeDistBounds,ampRatioLimit,findTrueSigma,iROI)

numPSFLocs = size(PSFLocs,1);
boxLength = 2*boxRadius+1;
PSFfits = NaN(numPSFLocs,18);

TxIdx = bsxfun(@plus,(-boxRadius:boxRadius)',PSFLocs(:,1)');
TyIdx = bsxfun(@plus,(-boxRadius:boxRadius)',PSFLocs(:,2)');

fitParam = zeros(numPSFLocs,9);

TxIdx(:,min(TxIdx)<1) = bsxfun(@plus,TxIdx(:,min(TxIdx)<1),(1-min(TxIdx(:,min(TxIdx)<1))));
TyIdx(:,min(TyIdx)<1) = bsxfun(@plus,TyIdx(:,min(TyIdx)<1),(1-min(TyIdx(:,min(TyIdx)<1))));
TxIdx(:,max(TxIdx)>cropWidth) = bsxfun(@minus,TxIdx(:,max(TxIdx)>cropWidth),max(TxIdx(:,max(TxIdx)>cropWidth)-cropWidth));
TyIdx(:,max(TyIdx)>cropHeight) = bsxfun(@minus,TyIdx(:,max(TyIdx)>cropHeight),max(TyIdx(:,max(TyIdx)>cropHeight)-cropHeight));

fitParam(:,3:8) = [PSFLocs(:,1) + templateLocs(PSFLocs(:,3),1)-(templateSize/2+0.5),PSFLocs(:,2) + templateLocs(PSFLocs(:,3),2)-(templateSize/2+0.5) ...
    PSFLocs(:,1) + templateLocs(PSFLocs(:,3),3)-(templateSize/2+0.5), PSFLocs(:,2) + templateLocs(PSFLocs(:,3),4)-(templateSize/2+0.5),repmat( mean(sigmaBounds),size(PSFLocs,1),2)];
withinBounds = round(fitParam(:,3)) > 0 & round(fitParam(:,4)) > 0 & round(fitParam(:,5)) > 0 & round(fitParam(:,6)) > 0 & round(fitParam(:,3)) < cropWidth & round(fitParam(:,5)) < cropWidth & round(fitParam(:,4)) < cropHeight & round(fitParam(:,6)) < cropHeight;
fitParam(withinBounds,1:2) = [data(sub2ind(size(data),round(fitParam(withinBounds,4)),round(fitParam(withinBounds,3)))),data(sub2ind(size(data),round(fitParam(withinBounds,6)),round(fitParam(withinBounds,5))))];


PSFfits(~withinBounds,1:8) = fitParam(~withinBounds,1:8);
PSFfits(~withinBounds,11) = -1000;

lowerBound = double([zeros(numPSFLocs,2), max([min(TxIdx);fitParam(:,3)'-4])', max([min(TyIdx);fitParam(:,4)'-4])', max([min(TxIdx);fitParam(:,5)'-4])', max([min(TyIdx);fitParam(:,6)'-4])', repmat(sigmaBounds(1),numPSFLocs,2), repmat(-median(median(bkgndImg)),numPSFLocs,1)]);
upperBound = double([zeros(numPSFLocs,2), min([max(TxIdx);fitParam(:,3)'+4])', min([max(TyIdx);fitParam(:,4)'+4])', min([max(TxIdx);fitParam(:,5)'+4])', min([max(TyIdx);fitParam(:,6)'+4])', repmat(sigmaBounds(2),numPSFLocs,2), repmat(median(median(bkgndImg)),numPSFLocs,1)]);

switch(fittingMethod)
    case 'LSQ with DG model'
        
    case 'MLE with DG model'
        %data is now in units of lambda (see Bewersdorf paper from
        %nature and methods)
        data2 = data + datavar;
        bkgndImg = double(bkgndImg);
        data2 = data2+bkgndImg;
        %bkgnd is now in photons and not counts!!
end

for b=1:numPSFLocs
    
    if ~all(all(FOVmask1(TyIdx(:,b),TxIdx(:,b))))
        PSFfits(b,1:8) = fitParam(b,1:8);
        PSFfits(b,11) = -1000;
    end
    if PSFfits(b,11) == -1000
        continue
    end
    
    upperBound(b,1:2) = repmat(1.2.*max(max(data(TyIdx(:,b),TxIdx(:,b)))),1,2);
%     upperBound(b,1:2) = repmat(max(max(data(TyIdx(:,b),TxIdx(:,b)))),1,2);
    xIdx = repmat(TxIdx(:,b),1,boxLength)';
    yIdx = repmat(TyIdx(:,b),1,boxLength);
    
    switch(fittingMethod)
        case 'LSQ with DG model'
            
            options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
            [PSFfits(b,1:8),~,residual,exitflag] = lsqnonlin(@(PARAM) ...
                f_doubleGaussianVector(PARAM,data(TyIdx(:,b),TxIdx(:,b)),0,xIdx,yIdx),...
                fitParam(b,1:8),lowerBound(b,1:8),upperBound(b,1:8),options);
            
            if findTrueSigma
                [PSFfits(b,7:8),~,~,~] = lsqnonlin(@(x)f_doubleGaussianVector([PSFfits(1:6) x],data(TyIdx(:,b),TxIdx(:,b)),0,xIdx,yIdx,fitParam(7:8),repmat(sigSearchBounds(1),1,2),repmat(sigSearchBounds(2),1,2),sigOptions));
            end
            %%divide by gain, convert to photons, for both
            fittedBkgndMean = mean2(bkgndImg(TyIdx(:,b),TxIdx(:,b)));
            PSFfits(b,9:11) = [fittedBkgndMean sum(abs(residual)) exitflag];
            
        case 'MLE with DG model'

            hMLE= @(PARAM)Likelihood(PARAM, data2(TyIdx(:,b),TxIdx(:,b)), bkgndImg(TyIdx(:,b),TxIdx(:,b)), datavar(TyIdx(:,b),TxIdx(:,b)),xIdx,yIdx);
            
            options = optimset('Display','off');
            [PSFfits(b,1:9), ~, exitflag] = fmincon(hMLE, fitParam(b,1:9), [], [], [], [], lowerBound(b,1:9), upperBound(b,1:9),  [], options);
            
            reconstructed = PSFfits(b,1).*exp( -((xIdx-PSFfits(b,3)).^2+(yIdx-PSFfits(b,4)).^2.) / (2.*PSFfits(b,7).^2))+PSFfits(b,2).*exp( -((xIdx-PSFfits(b,5)).^2+(yIdx-PSFfits(b,6)).^2.) / (2.*PSFfits(b,8).^2))+PSFfits(b,9);
%             residual = sum(sum(abs(reconstructed+bkgndImg(TyIdx(:,b),TxIdx(:,b))-data(TyIdx(:,b),TxIdx(:,b)))));
            residual = sum(sum(abs(reconstructed-data(TyIdx(:,b),TxIdx(:,b))))); %attempting to fix fiterror. JR 11/30/2016

            
            PSFfits(b,10:11) = [residual exitflag];
    end
    PSFfits(b,15) = sum(sum(data(TyIdx(:,b),TxIdx(:,b))));
end
% shift coordinates relative to entire dataset (not just ROI)
PSFfits(:,3) = PSFfits(:,3) + ROI(1)-1+iROI(1)-1;
PSFfits(:,4) = PSFfits(:,4) + ROI(2)-1+iROI(2)-1;
PSFfits(:,5) = PSFfits(:,5) + ROI(1)-1+iROI(1)-1;
PSFfits(:,6) = PSFfits(:,6) + ROI(2)-1+iROI(2)-1;
% Calculate midpoint between two Gaussian spots
% shift coordinates relative to entire dataset (not just ROI) and
% convert from pixels to nm
PSFfits(:,12) = (PSFfits(:,3)+PSFfits(:,5))/2*nmPerPixel;
PSFfits(:,13) = (PSFfits(:,4)+PSFfits(:,6))/2*nmPerPixel;

% Below is the calculation of the angle of the two lobes.
% Remember that two vertical lobes is focal plane because camera
% outputs data that is rotated. Therefore, we want y2>y1 for all
% angle calculations (so that -90<=angle<=90, and we use swap
% the use of x and y for the atan2 calculation.
%     x1 = fitParam(3);
%     x2 = fitParam(5);
%     y1 = fitParam(4);
%     y2 = fitParam(6);
%     % swap if y1>y2
%     if (y1 > y2)
%         tx = x1; ty = y1;
%         x1 = x2; y1 = y2;
%         x2 = tx; y2 = ty;
%         %clear tx ty;
%     end
%Finds the angle
flop = PSFfits(:,6)>PSFfits(:,4);
PSFfits(flop,14) = atan2(-(PSFfits(flop,5)-PSFfits(flop,3)),PSFfits(flop,6)-PSFfits(flop,4)) * 180/pi;
PSFfits(~flop,14) = atan2(-(PSFfits(~flop,3)-PSFfits(~flop,5)),PSFfits(~flop,4)-PSFfits(~flop,6)) * 180/pi;

%Below is a way to count the photons. It integrates the box and
%subtracts the boxarea*offset from the fit. It is inherently flawed
%if there happens to be bright pixels inside of the fitting
%region

%     totalCounts = sum(sum(data(yIdx(:,1),xIdx(1,:))));
%     totalCounts = zeros(numPSFLocs,1)
% % for b = 1:numPSFLocs
% %     switch(fittingMethod)
% %         case 'LSQ with DG model'
% %             PSFfits(b,15) = sum(sum(data(TyIdx(:,b),TxIdx(:,b))));
% %         case 'MLE with DG model'
% %             %%We do not want to include bkgnd photons in the
% %             %%total count, so we subtract. In the case of LSQ
% %             %%they were already taken out of data
% %             PSFfits(b,15) = sum(sum(data(TyIdx(:,b),TxIdx(:,b))))- sum(sum(bkgndImg(TyIdx(:,b),TxIdx(:,b))))-sum(sum(datavar(TyIdx(:,b),TxIdx(:,b))));
% %     end
% % end
%  PSFfits(b,15) = (totalCounts-(2*boxRadius+1)^2*bkgndMean)*conversionFactor;  % bkgndMean = 0
%Do not multiply by conversion gain: already in terms of
%photons if LSQ and lambda if MLE
%     PSFfits(:,15) = totalCounts';  % bkgndMean = 0
% PSFfits(b,15) = totalCounts*conversionFactor;

%The interlobe distance
PSFfits(:,16) = sqrt((PSFfits(:,3)-PSFfits(:,5)).^2 + (PSFfits(:,4)-PSFfits(:,6)).^2);
%Amplitude Ratio
PSFfits(:,17) = abs(PSFfits(:,1) - PSFfits(:,2))./sum(PSFfits(:,1:2),2);
% Gaussian width Ratio
PSFfits(:,18) = abs(PSFfits(:,7) - PSFfits(:,8))./sum(PSFfits(:,7:8),2);

%% Now evaluate the goodness of the fits

% Conditions for fits (play with these):
% (1) Amplitude of both lobes > 0
% (2) All locations x1,y1, x2,y2 lie inside area of small box
% (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
% (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
% (5) Make sure amplitudes are within 100% of one another
% (6) Make sure totalFitError/(total number of photons) < 1.05

% absolute amplitude > 0?
PSFfits(PSFfits(:,1)<0,11) = -1001;
PSFfits(PSFfits(:,2)<0,11) = -1001;
% peaks inside box?
PSFfits(PSFfits(:,3)<(TxIdx(1,:)+ROI(1)-1+iROI(1)-1)',11) = -1002;
PSFfits(PSFfits(:,4)<(TyIdx(1,:)+ROI(2)-1+iROI(2)-1)',11) = -1002;
PSFfits(PSFfits(:,5)<(TxIdx(1,:)+ROI(1)-1+iROI(1)-1)',11) = -1002;
PSFfits(PSFfits(:,6)<(TyIdx(1,:)+ROI(2)-1+iROI(2)-1)',11) = -1002;
PSFfits(PSFfits(:,3)>(TxIdx(boxLength,:)+ROI(1)-1+iROI(1)-1)',11) = -1002;
PSFfits(PSFfits(:,4)>(TyIdx(boxLength,:)+ROI(2)-1+iROI(2)-1)',11) = -1002;
PSFfits(PSFfits(:,5)>(TxIdx(boxLength,:)+ROI(1)-1+iROI(1)-1)',11) = -1002;
PSFfits(PSFfits(:,6)>(TyIdx(boxLength,:)+ROI(2)-1+iROI(2)-1)',11) = -1002;
% absolute sigma size for either lobe within bounds?
if findTrueSigma
    PSFfits(PSFfits(:,7)-0.001 <= sigSearchBounds(1),11) = -1003;
    PSFfits(PSFfits(:,8)-0.001 <= sigSearchBounds(1),11) = -1003;
    PSFfits(PSFfits(:,7)+0.001 >= sigSearchBounds(2),11) = -1003;
    PSFfits(PSFfits(:,8)+0.001 >= sigSearchBounds(2),11) = -1003;
else
    PSFfits(PSFfits(:,7) <= sigSearchBounds(1),11) = -1003;
    PSFfits(PSFfits(:,8) <= sigSearchBounds(1),11) = -1003;
    PSFfits(PSFfits(:,7) >= sigSearchBounds(2),11) = -1003;
    PSFfits(PSFfits(:,8) >= sigSearchBounds(2),11) = -1003;
end
% sigma ratio of lobes less than limit?
PSFfits(PSFfits(:,18) > sigmaRatioLimit,11) = -1004;
% interlobe distance within bounds?
PSFfits(PSFfits(:,16) > lobeDistBounds(2),11) = -1005;
PSFfits(PSFfits(:,16) < lobeDistBounds(1),11) = -1005;
% amplitude ratio of lobes less than limit?
PSFfits(PSFfits(:,17) > ampRatioLimit,11) = -1006;

end