function refit = RefitErrors_2(refit,TxIdx,TyIdx,data,bkgndImg,datavar,fittingMethod,numBeads,nmPerPixel,horizontal,boxLength,a,c, sigmaBounds,lB,uB,sigmaRatioLimit,lobeDistBounds,ampRatioLimit, cropHeight, cropWidth)
ibadFits = 0;
mBeads = nan(size(refit)); sBeads = nan(size(refit));
reProcess = refit(:,13)<0;
for b = 1:numBeads
    [~,iBeads] = sort((refit(:,14)-refit(b,14)).^2 + (refit(:,15)-refit(b,15)).^2);
    rBeads = iBeads(~reProcess(iBeads));
    rBeads(rBeads== b) = [];
    gBeads = refit(rBeads(1:min(20,length(rBeads))),:);
    gBeads(:,5:8) = gBeads(:,5:8) - repmat([TxIdx(1,rBeads(1:min(20,length(rBeads))))',TyIdx(1,rBeads(1:min(20,length(rBeads))))'],1,2);
    mBeads(b,:) = mean(gBeads);
    sBeads(b,:) =  std(gBeads);
end

aBeads = refit(:,13) < -1000 |...
    refit(:,3) > mBeads(:,3)+3*sBeads(:,3)|...
    refit(:,4) > mBeads(:,4)+3*sBeads(:,4)|...
    refit(:,3) < mBeads(:,3)-3*sBeads(:,3)|...
    refit(:,4) < mBeads(:,4)-3*sBeads(:,4)|...
    refit(:,5) - TxIdx(1,:)' > mBeads(:,5)+3*sBeads(:,5)|...
    refit(:,6) - TyIdx(1,:)' > mBeads(:,6)+3*sBeads(:,6)|...
    refit(:,7) - TxIdx(1,:)' > mBeads(:,7)+3*sBeads(:,7)|...
    refit(:,8) - TyIdx(1,:)' > mBeads(:,8)+3*sBeads(:,8)|...
    refit(:,5) - TxIdx(1,:)' < mBeads(:,5)-3*sBeads(:,5)|...
    refit(:,6) - TyIdx(1,:)' < mBeads(:,6)-3*sBeads(:,6)|...
    refit(:,7) - TxIdx(1,:)' < mBeads(:,7)-3*sBeads(:,7)|...
    refit(:,8) - TyIdx(1,:)' < mBeads(:,8)-3*sBeads(:,8)|...
    refit(:,16) >  mBeads(:,16)+3*sBeads(:,16) |...
    refit(:,16) <  mBeads(:,16)-3*sBeads(:,16) |...
    refit(:,19) >  mBeads(:,19)+3*sBeads(:,19) |...
    refit(:,19) <  mBeads(:,19)-3*sBeads(:,19);
bBeads = refit(aBeads & refit(:,13) ~= -1000,:);
badFits = size(bBeads,1);
while badFits > ibadFits
    badFits = size(bBeads,1);
    cBeads = mBeads(bBeads(:,2),3:8) + [zeros(badFits,2),repmat([TxIdx(1,bBeads(:,2))',TyIdx(1,bBeads(:,2))'],1,2)];
    refitS = cell(badFits,1);
    if badFits > 24
        parfor d = 1:size(bBeads,1)
            b = bBeads(d,2);
            
            xIdx = repmat(TxIdx(:,b),1,boxLength)';
            yIdx = repmat(TyIdx(:,b),1,boxLength);
            fitParam(1:6) = cBeads(d,:);
            fitParam(7) = mean(sigmaBounds);
            fitParam(8) = mean(sigmaBounds);
            lowerBound = lB(b,1:8);
            upperBound = uB(b,1:8);
            upperBound(1:2) = repmat(1.2*max(max(data(yIdx(:,1),xIdx(1,:)))),1,2);
            switch(fittingMethod)
                case 'LSQ with DG model'
                    %% Fit with lsqnonlin
                    options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on','Display','off');
                    [fitParam,~,residual,exitflag] = lsqnonlin(@(x) ...
                        f_doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),0,xIdx,yIdx),...
                        fitParam,lowerBound,upperBound,options);
                    %         bkgndCnts = bkgndImg(round((fitParam(4)+fitParam(6))/2),round((fitParam(3)+fitParam(5))/2));
                    bkgndCnts = 0;
                    %%Note we need 9 params for MLE, most likely just get rid of
                    %%bkgnd counts, no residual,no exit flag(just set as 1?)
                    refitS{d,1} = [a+c b fitParam bkgndCnts sum(abs(residual)) exitflag];
                case 'MLE with DG model'
                    lowerBound(9) = -median(median(bkgndImg));
                    upperBound(9) = median(median(bkgndImg));
                    fitParam(9) = 0;
                    hMLE= @(fitParam)Likelihood( fitParam, data(yIdx(:,1),xIdx(1,:)), bkgndImg(yIdx(:,1),xIdx(1,:)), datavar(yIdx(:,1),xIdx(1,:)), xIdx,yIdx);
                    
                    options = optimset('Display','off');
                    [fitParam, ~, exitflag] = fmincon(hMLE, fitParam, [], [], [], [], lowerBound, upperBound, [], options);
                    
                    %%Note we need 9 params for MLE, most likely just get rid of
                    %%bkgnd counts, no residual,no exit flag(just set as 1?)
                    %%Changed to reflect MLE
                    %%could use fval instead of NaN
                    refitS{d,1} = [a+c b fitParam nan exitflag];
            end
        end
    else
        for d = 1:size(bBeads,1)
            b = bBeads(d,2);
            
            xIdx = repmat(TxIdx(:,b),1,boxLength)';
            yIdx = repmat(TyIdx(:,b),1,boxLength);
            fitParam(1:6) = cBeads(d,:);
            fitParam(7) = mean(sigmaBounds);
            fitParam(8) = mean(sigmaBounds);
            lowerBound = lB(b,1:8);
            upperBound = uB(b,1:8);
            upperBound(1:2) = repmat(1.2*max(max(data(yIdx(:,1),xIdx(1,:)))),1,2);
            switch(fittingMethod)
                case 'LSQ with DG model'
                    %% Fit with lsqnonlin
                    options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on','Display','off');
                    [fitParam,~,residual,exitflag] = lsqnonlin(@(x) ...
                        f_doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),0,xIdx,yIdx),...
                        fitParam,lowerBound,upperBound,options);
                    %         bkgndCnts = bkgndImg(round((fitParam(4)+fitParam(6))/2),round((fitParam(3)+fitParam(5))/2));
                    bkgndCnts = 0;
                    %%Note we need 9 params for MLE, most likely just get rid of
                    %%bkgnd counts, no residual,no exit flag(just set as 1?)
                    refitS{d,1} = [a+c b fitParam bkgndCnts sum(abs(residual)) exitflag];
                case 'MLE with DG model'
                    lowerBound(9) = -median(median(bkgndImg));
                    upperBound(9) = median(median(bkgndImg));
                    fitParam(9) = 0;
                    hMLE= @(fitParam)Likelihood( fitParam, data(yIdx(:,1),xIdx(1,:)), bkgndImg(yIdx(:,1),xIdx(1,:)), datavar(yIdx(:,1),xIdx(1,:)), xIdx,yIdx);
                    
                    options = optimset('Display','off');
                    [fitParam, ~, exitflag] = fmincon(hMLE, fitParam, [], [], [], [], lowerBound, upperBound, [], options);
                    
                    %%Note we need 9 params for MLE, most likely just get rid of
                    %%bkgnd counts, no residual,no exit flag(just set as 1?)
                    %%Changed to reflect MLE
                    %%could use fval instead of NaN
                    refitS{d,1} = [a+c b fitParam nan exitflag];
            end
        end
    end
    refit(bBeads(:,2),1:13) = cell2mat(refitS);
    for b = 1:numBeads
        if min(TxIdx(:,b)) < 1 || max(TyIdx(:,b)) > cropHeight || max(TxIdx(:,b)) > cropWidth || min(TyIdx(:,b)) < 1
            refit(b,13) = -1000;
            refit(b,17) = NaN;
            continue
        end
        %Below is a way to count the photons. It integrates the box and
        %subtracts the boxarea*offset from the fit. It is inherently flawed
        %if there happen to be bright pixels inside of the fitting region.
        %Units of photons for LSQ nonlin
        refit(b,17) = sum(sum(data(TyIdx(:,b),TxIdx(:,b))));
    end

    refit(:,14) = (refit(:,5)+refit(:,7))/2*nmPerPixel;
    refit(:,15) = (refit(:,6)+refit(:,8))/2*nmPerPixel;
    if ~horizontal
        flop = refit(:,8)>refit(:,6);
        refit(flop,16) = atan2(-(refit(flop,7)-refit(flop,5)),refit(flop,8)-refit(flop,6)) * 180/pi;
        refit(~flop,16) = atan2(-(refit(~flop,5)-refit(~flop,7)),refit(~flop,6)-refit(~flop,8)) * 180/pi;
    elseif horizontal
        flop = refit(:,7)>refit(:,5);
        refit(flop,16) = atan2(-(refit(flop,8)-refit(flop,6)),refit(flop,7)-refit(flop,5)) * 180/pi;
        refit(~flop,16) = atan2(-(refit(~flop,6)-refit(~flop,8)),refit(~flop,5)-refit(~flop,7)) * 180/pi;
    end
    
    %The interlobe distance
    refit(:,19) = sqrt((refit(:,5)-refit(:,7)).^2 + (refit(:,6)-refit(:,8)).^2);
    %Amplitude Ratio
    refit(:,20) = abs(refit(:,3) - refit(:,4))./sum(refit(:,3:4),2);
    % Gaussian width Ratio
    refit(:,21) = abs(refit(:,9) - refit(:,10))./sum(refit(:,9:10),2);
    
    %% Now evaluate the fits
    % Conditions for fits (play with these):
    % (1) Amplitude of both lobes > 0
    % (2) All locations x1,y1, x2,y2 lie inside area of small box
    % (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
    % (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
    % (5) Make sure amplitudes are within 100% of one another
    % (6) Make sure totalFitError/(total number of photons) < 1.05 (not
    %     enabled at the present time)
    
    
    % absolute amplitude > 0?
    refit(refit(:,3)<0,13) = -1001;
    refit(refit(:,4)<0,13) = -1001;
    % peaks inside box?
    % sigma ratio of lobes less than limit?
    refit(refit(:,21) > sigmaRatioLimit,13) = -1004;
    % interlobe distance within bounds?
    refit(refit(:,19) > lobeDistBounds(2),13) = -1005;
    refit(refit(:,19) < lobeDistBounds(1),13) = -1005;
    % amplitude ratio of lobes less than limit?
    refit(refit(:,20) > ampRatioLimit,13) = -1006;

    mBeads = nan(size(refit)); sBeads = nan(size(refit));
    reProcess = refit(:,13)<0;
    for b = 1:numBeads
        [dBeads,iBeads] = sort((refit(:,14)-refit(b,14)).^2 + (refit(:,15)-refit(b,15)).^2);
        rBeads = iBeads(~reProcess(iBeads));
        rBeads(rBeads== b) = [];
        gBeads = refit(rBeads(1:min(20,length(rBeads))),:);
        gBeads(:,5:8) = gBeads(:,5:8) - repmat([TxIdx(1,rBeads(1:min(20,length(rBeads))))',TyIdx(1,rBeads(1:min(20,length(rBeads))))'],1,2);
        mBeads(b,:) = mean(gBeads);
        sBeads(b,:) =  std(gBeads);
    end
    
    aBeads = refit(:,13) < -1000 |...
        refit(:,3) > mBeads(:,3)+3*sBeads(:,3)|...
        refit(:,4) > mBeads(:,4)+3*sBeads(:,4)|...
        refit(:,3) < mBeads(:,3)-3*sBeads(:,3)|...
        refit(:,4) < mBeads(:,4)-3*sBeads(:,4)|...
        refit(:,5) - TxIdx(1,:)' > mBeads(:,5)+3*sBeads(:,5)|...
        refit(:,6) - TyIdx(1,:)' > mBeads(:,6)+3*sBeads(:,6)|...
        refit(:,7) - TxIdx(1,:)' > mBeads(:,7)+3*sBeads(:,7)|...
        refit(:,8) - TyIdx(1,:)' > mBeads(:,8)+3*sBeads(:,8)|...
        refit(:,5) - TxIdx(1,:)' < mBeads(:,5)-3*sBeads(:,5)|...
        refit(:,6) - TyIdx(1,:)' < mBeads(:,6)-3*sBeads(:,6)|...
        refit(:,7) - TxIdx(1,:)' < mBeads(:,7)-3*sBeads(:,7)|...
        refit(:,8) - TyIdx(1,:)' < mBeads(:,8)-3*sBeads(:,8)|...
        refit(:,16) >  mBeads(:,16)+3*sBeads(:,16) |...
        refit(:,16) <  mBeads(:,16)-3*sBeads(:,16) |...
        refit(:,19) >  mBeads(:,19)+3*sBeads(:,19) |...
        refit(:,19) <  mBeads(:,19)-3*sBeads(:,19);
    bBeads = refit(aBeads & refit(:,13) ~= -1000,:);
    ibadFits = size(bBeads,1);
end
end