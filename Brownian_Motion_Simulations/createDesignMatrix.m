%createDesignMatrix
%This function creates a matrix of curves from the scatteredInterpolant F,
%for linear fitting of curves

[dataFile3, dataPath3] = uigetfile({'*.mat';'*.*'},'Open interpolated curves used for fitting','MultiSelect', 'on');
load([dataPath3 dataFile3]);

xq = [0:0.01:20] ;    % vector of query points, needs finer resolution at low d



resolution = 0.05;
dMax = 15; % max diffusion coefficient
dRange = 1:resolution:dMax;


C = zeros(length(xq),length(dRange));

figure;
for i = 1:length(dRange)
    yq = ones(1,length(xq))*dRange(i);
    C(:,i) = F(xq,yq);
%     C(:,i) = F(xq,yq)./sum(F(xq,yq));

    plot(xq,C(:,i));
    hold on;
end
    save(['designMatrix_' num2str(resolution) '_20.mat'], 'C','xq', 'resolution', 'dMax', 'dRange')

z = 1;