% Julian Rocha

%addDiffusionCoefficients
% This function inputs several separate files with diffusion coefficient
% information and combines them. This is typically used to combine data for
% several field's-of-view (FOV). 

[dataFile, dataPath] = uigetfile({'*.mat';'*.*'},'Open all files with diffusion coefficients','MultiSelect', 'on');

totalDiffusion = [];
numCellsTotal = 0;
if iscell(dataFile)
    for a = 1:length(dataFile)
        load ([dataPath,dataFile{a}]);
        totalDiffusion = vertcat(totalDiffusion,diffusionCoefficients);
        numCellsTotal = numCellsTotal + numCells;
    end
else
    load ([dataPath,dataFile]);
    totalDiffusion = vertcat(totalDiffusion,diffusionCoefficients);
    numCellsTotal = numCellsTotal + numCells;
end

save('totalDiffusion.mat');

% %plot histogram of results
% figure;
% hold on;
% title(['Diffusion Coefficients: ' num2str(length(totalDiffusion)) ' Tracks, ' num2str(numCellsTotal) ' Cells ']);
% xlabel('Diffusion Coefficient (µm^2/s)');
% ylabel('Probability Density');
% edges = 0:0.1:10;
% binSpacing = edges(2)-edges(1);
% histogram_counts = histcounts(totalDiffusion,(edges+binSpacing/2));
% pdf_distribution = histogram_counts./(sum(histogram_counts));
% xlim([0 max(edges)]);
% bar1 = bar(edges(2:end),pdf_distribution); 