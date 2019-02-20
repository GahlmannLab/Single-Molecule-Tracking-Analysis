function [F] = Bspline(order, numPoints, xData, yData)

% function spline(n,order)
%
% Plots the B-slpine-curve of n control-points.
% The control points can be chosen by clicking
% with the mouse on the figure.
%
% COMMAND:  spline(n,order)
% INPUT:    n     Number of Control-Points
%           order Order of B-Splines
%                 Argnument is arbitrary
%                 default: order = 4

%% B spline a 2D data
% [dataFile, dataPath] = uigetfile({'*.mat';'*.*'},'MultiSelect','on','Open .mat calibration file of NHA');
% 
% if ischar(dataFile)==1
%      dataFile = {dataFile};
% end
% numFiles = length(dataFile);
% for f = 1:numFiles;
% load ([dataPath,dataFile{f}]);
% end

if size(xData,1) == 1
    xData = xData';
end
if size(yData,1) == 1
    yData = yData';
end

xy = [xData,yData];
n = length(xData);
% order=20;
BsplineFit = [];
 
% nplot = 100;
% t = linspace(0,1,nplot);

if (n < order)
	display([' !!! Error: Choose n >= order=',num2str(order),' !!!']);
	return;
end

%% Plot curves for x shift
% for j = 1:m
%     xy = [angles_NHA_good,meanX_good];
    T = linspace(0,1,n-order+2);
    y = linspace(0,1,numPoints);
    BsplineFit = DEBOOR(T,xy,y,order);
% end

% %% plot results
% figure('Position',([100 100 1000 700]));
% hold on; box on;
% set(gca,'Fontsize',16);
% % for k = 1:m
% %     h_x = plot(p_spl_x{k}(:,1),p_spl_x{k}(:,2),'r-','LineWidth',0.5,'color',rand(1,3));
% % end
% plot(BsplineFit(:,1),BsplineFit(:,2),'r-','LineWidth',0.5); % plot results after spline
% plot(xy(:,1),xy(:,2),'bo','LineWidth',0.5); % plot original points
% 
% % xlim([-150 100]);
% xlabel('x');
% ylabel('y');
% title('B Spline Fit');

%% Create and Gridded Interpolant that can be queried at any position
%JR 11/1/2018. Added so that the gridded interpolant can interpolate
%correctly near the upper limit of the range
padSize = ceil(xData(end)-max(BsplineFit(:,1)))-1;
padStart = ceil(max(BsplineFit(:,1)))+1;
if padSize <= 0
    padSize = 1;
    padStart = xData(end);
end

tempX = BsplineFit(:,1);
tempY = BsplineFit(:,2);
tempY = padarray(tempY,padSize, 1,'post');
for m = 1:padSize
    tempX = padarray(tempX,1,padStart,'post');
    padStart = padStart + 1;
end
BsplineFit = [tempX, tempY];
% BsplineFit = [BsplineFit; max(xData),1];


F = griddedInterpolant(BsplineFit(:,1),BsplineFit(:,2), 'spline');
% plot(xData,F(xData),'ro','LineWidth',0.5); % plot original points


end