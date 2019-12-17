% Julian Rocha

interpxVals = [];
yVals = [];
vVals = []; 


%range of the distribution
xq = 0:0.001:25;

%diffusion coefficients being generated
diffVals = [1:1:15];

load('allDiffCoeff_1-15.mat');

allVals = allDiffCoeff;


alpha = 0.300;
%% Plot all eCDFs
curvesTemp = cell(length(diffVals),1);
curvesTempX = cell(length(diffVals),1);
curvesInterp = cell(length(diffVals),1);
Curves = zeros(length(diffVals),length(xq));

for i = 1: length(allVals)

    %calculated ecdf
    ecdf(allVals{i})
    [curvesTemp{i}, curvesTempX{i}] =  ecdf(allVals{i});
    
    %pad curvesTemp on both sides so that the Bspline interpolant can be
    %queried outside of the range of the ecdf curve
    padSize1 = floor(min(curvesTempX{i})/0.001);
    padStart1 = padSize1*0.001;
    
    curvesTemp{i} = padarray(curvesTemp{i},padSize1,0,'pre');
%     curvesTempX{i} = padarray(curvesTempX{i},padSize1,0,'pre');
    for j = 1:padSize1
    curvesTempX{i} = padarray(curvesTempX{i},1,padStart1,'pre');
    padStart1 = padStart1 - 0.001;
    end
%     curvesTemp{i} = padarray(curvesTemp{i},1,0,'pre');
%     curvesTempX{i} = padarray(curvesTempX{i},1,0,'pre');
    
    
    padSize2 = ceil(xq(end)-max(curvesTempX{i}));
%     padSize2 = ceil(25-max(curvesTempX{i}));
%     padSize2 = ceil(30-max(curvesTempX{i}));
    padStart2 = ceil(max(curvesTempX{i}));
    if padSize2 < 0 
        padSize2 = 1;
        padStart2 = xq(end);
%         padStart2 = 25;
%         padStart2 = 30;
    end
   
    
    curvesTemp{i} = padarray(curvesTemp{i},padSize2,1,'post');
    for j = 1:padSize2
    curvesTempX{i} = padarray(curvesTempX{i},1,padStart2,'post');
    padStart2 = padStart2 + 1;
    end
    %end pad curvesTemp
    
        %interpolate ecdf so that all curves have the same number of points,
    %necessary for scatteredInterpolant function later
%     Curves(i,:) = interp1(curvesTempX{i},curvesTemp{i},xq,'spline');
    order = 10;
%     order = 40;
%     numPoints = length(curvesTempX{i});
    numPoints = 100000;
    curvesInterp{i} = Bspline(order,numPoints,curvesTempX{i},curvesTemp{i});
    
    Curves(i,:) = curvesInterp{i}(xq);
%     %correct Curves for out of range values
%     tempCurve = Curves(i,:);
%     tempCurve(xq < min(curvesTempX{i})) = 0;
%     tempCurve(xq > max(curvesTempX{i})) = 1;
%     Curves(i,:) = tempCurve; 
    i
        
end
% hold off
close all

figure; 
for i = 1:length(allVals)
    hold on  
    plot(xq,Curves(i,:))
end
hold off


xVals = repmat(xq,length(diffVals),1);
xVals = reshape(xVals,[],1);
yVals = repmat(diffVals,length(xq),1);
yVals = yVals';
yVals = reshape(yVals,[],1);
vVals = reshape(Curves,[],1);
F = scatteredInterpolant(xVals,yVals,vVals, 'natural');

z = 1;
%save the variable F, use for fitting


%%
% queryCurves = 0.01:0.01:20;
% % [xq1, yq1] = meshgrid(xq,queryCurves);
% tic
% [xq1, yq1] = meshgrid(0.01:0.01:2,2:0.001:4);
% vq1 = F(xq1,yq1);
% 
% figure;
% imagesc(vq1, [0 1]);
% toc

