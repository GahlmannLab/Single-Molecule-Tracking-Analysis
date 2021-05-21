function PSFfitsS = PSFfitting_calDHPSF(data,fittingMethod,cropHeight,cropWidth,TxIdx,TyIdx,sigmaBounds,moleLocs,dataBlur,datavar,bkgndImg,boxLength,a,c,b,uB,lB)

% create indices to use for fitting
xIdx = repmat(TxIdx(:,b),1,boxLength)';
yIdx = repmat(TyIdx(:,b),1,boxLength);
% make sure indices are inside image

% find two largest peaks in box
[tempY, tempX] = ind2sub([boxLength boxLength], ...
    find(imregionalmax(dataBlur(yIdx(:,1),xIdx(1,:)))));
tempX = tempX + min(xIdx(:))-1;
tempY = tempY + min(yIdx(:))-1;

temp = sortrows([tempX tempY data(sub2ind([cropHeight cropWidth], ...
    tempY,tempX))],-3);


% set initial fitting parameters
% [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
% if two peaks aren't found, then use previous fitting parameters
switch(fittingMethod)
    case 'LSQ with DG model'
        if size(temp,1)>=2
            fitParam(3) = temp(1,1);
            fitParam(4) = temp(1,2);
            fitParam(5) = temp(2,1);
            fitParam(6) = temp(2,2);
            fitParam(1) = temp(1,3);
            fitParam(2) = temp(2,3);
            fitParam(7) = mean(sigmaBounds);
            fitParam(8) = mean(sigmaBounds);
        end
        lowerBound = lB(b,1:8);
        upperBound = uB(b,1:8);
        upperBound(1:2) = repmat(1.2*max(max(data(yIdx(:,1),xIdx(1,:)))),1,2);
        
        %% Fit with lsqnonlin
        options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on','Display','off');
        [fitParam,~,residual,exitflag] = lsqnonlin(@(x) ...
            f_doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),0,xIdx,yIdx),...
            fitParam,lowerBound,upperBound,options);
        %         bkgndCnts = bkgndImg(round((fitParam(4)+fitParam(6))/2),round((fitParam(3)+fitParam(5))/2));
        bkgndCnts = 0;
        %%Note we need 9 params for MLE, most likely just get rid of
        %%bkgnd counts, no residual,no exit flag(just set as 1?)
        PSFfitsS = [a+c b fitParam bkgndCnts sum(abs(residual)) exitflag];
    case 'MLE with DG model'
        if size(temp,1)>=2
            fitParam(3) = temp(1,1);
            fitParam(4) = temp(1,2);
            fitParam(5) = temp(2,1);
            fitParam(6) = temp(2,2);
            fitParam(1) = temp(1,3);
            fitParam(2) = temp(2,3);
            fitParam(7) = mean(sigmaBounds);
            fitParam(8) = mean(sigmaBounds);
            fitParam(9) = 0;
        else
            fitParam(3) = moleLocs(b,1)+4;
            fitParam(4) = moleLocs(b,2)+4;
            fitParam(5) = moleLocs(b,1)-4;
            fitParam(6) = moleLocs(b,2)-4;
            fitParam(1) = data(fitParam(4),fitParam(3));
            fitParam(2) = data(fitParam(6),fitParam(5));
            fitParam(7) = mean(sigmaBounds);
            fitParam(8) = mean(sigmaBounds);
            fitParam(9) = 0;
        end
        %modified by Alecia Achimovich 5/21/2021 - If median value of
        %background image is negative, the result will be lower bound that
        %exceeds the upper bound (lower bound will be positive, and upper
        %bound will be negative. fmincon will error out.
        
        %lowerBound = [lB(b,1:8), -median(median(bkgndImg))]; 
        %upperBound = [uB(b,1:8), median(median(bkgndImg))];
        
        %Crop corners of FOV to determine median.
        ROImed = [data(1:10,1:10),data(cropHeight-9:cropHeight,1:10);data(1:10,cropWidth-9:cropWidth),data(cropHeight-9:cropHeight,cropWidth-9:cropWidth)];
        lowerBound = [lB(b,1:8), 0]; 
        upperBound = [uB(b,1:8),  1.5*median(median(ROImed))];
        upperBound(1:2) = repmat(1.2*max(max(data(yIdx(:,1),xIdx(1,:)))),1,2);
        %end modifications by Alecia Achimovich
        
        upperBound(1:2) = repmat(1.2*max(max(data(yIdx(:,1),xIdx(1,:)))),1,2);
        
        hMLE= @(fitParam)Likelihood( fitParam, data(yIdx(:,1),xIdx(1,:)), bkgndImg(yIdx(:,1),xIdx(1,:)), datavar(yIdx(:,1),xIdx(1,:)), xIdx,yIdx);
        
        options = optimset('Display','off');
        [fitParam, ~, exitflag] = fmincon(hMLE, fitParam, [], [], [], [], lowerBound, upperBound, [], options);
        
        %%Note we need 9 params for MLE, most likely just get rid of
        %%bkgnd counts, no residual,no exit flag(just set as 1?)
        %%Changed to reflect MLE
        %%could use fval instead of NaN
        PSFfitsS = [a+c b fitParam NaN exitflag];
end
