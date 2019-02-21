function meanAngles_NHA = LocalElimination2(absLocs,meanAngles,outputFilePrefix)
meanAngles_NHA = meanAngles;
numBeads = size(meanAngles_NHA,1);
steps = 1:size(meanAngles_NHA,2);
[~,iBeads] = sort(bsxfun(@hypot, bsxfun(@minus,absLocs(:,1),absLocs(:,1)'), bsxfun(@minus,absLocs(:,2),absLocs(:,2)')));
[~,BOrder] = sort(bsxfun(@hypot, bsxfun(@minus,absLocs(:,1),1024*108), bsxfun(@minus,absLocs(:,2),1024*108)));
BadBeads = max(sum(isnan(meanAngles_NHA(:))),1);
iBadBeads = 0;
% fig = figure('Position',[100,100,1600,800]); 
c = 1;
while iBadBeads < BadBeads
    iBadBeads = BadBeads;
    goodBeads = zeros(numBeads,1);
    goodBeads(1:(20+1),1) = BOrder(1:(20+1));
    for b = 1:numBeads
        a = BOrder(b);
        gBeads = iBeads(ismember(iBeads(:,a),goodBeads),a);
        gBeads(gBeads == a) = [];
%         gBeads(sum(isnan(meanAngles_NHA(gBeads,:)),2)>6) = [];
        bspace = 6.*nanmean(nanstd(meanAngles_NHA(gBeads(1:10),:)));
        mL = nanmedian(meanAngles_NHA(gBeads(1:10),:));
%         mL(isnan(mL)) = nanmedian(meanAngles_NHA(:,isnan(mL)))-nanmean(nanmedian(meanAngles_NHA)-meanAngles_NHA(a,:));
        mL = interp1(steps(~isnan(mL)),mL(~isnan(mL)),steps,'pchip',-90);
        uL = mL + bspace;
        lL = mL - bspace;
%         
%         if any(meanAngles_NHA(a,:)>uL | meanAngles_NHA(a,:)<lL)
%             subplot(1,2,1); hold off;
%             plot(nanmean(meanAngles_NHA(gBeads(1:10),:)),'color','black');hold on;
%             plot(meanAngles_NHA(gBeads(1:10),:)','--');
%             plot(uL,':','color','red'); plot(lL,':','color','red');
%             plot(meanAngles_NHA(a,:),'x','color','black');
%             plot(steps(meanAngles_NHA(a,:)>uL | meanAngles_NHA(a,:)<lL),meanAngles_NHA(a,meanAngles_NHA(a,:)>uL | meanAngles_NHA(a,:)<lL),'o','color','blue');
%             title(['Beads: ' num2str(a)]);
%             subplot(1,2,2); hold off;
%             plot(absLocs(gBeads,1),absLocs(gBeads,2),'.','color','blue');hold on;
%             plot(absLocs(gBeads(1:10),1),absLocs(gBeads(1:10),2),'o','color','red');
%             plot(absLocs(a,1),absLocs(a,2),'x','color','black');
%             xlim([0,108*2048]);ylim([0,108*2048]);
%             title(['Pass: ' num2str(c)]);
%             saveas(fig,[outputFilePrefix 'b' num2str(a) '_p' num2str(c) '_onesided_mean.png']);
%         end
        meanAngles_NHA(a,(meanAngles_NHA(a,:)>uL | meanAngles_NHA(a,:)<lL)) = nan;
        goodBeads(b,1) = a;
    end
    BadBeads = sum(isnan(meanAngles_NHA(:)));
    c = c+1;
end
end