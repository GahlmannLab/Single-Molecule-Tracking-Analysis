function [curveInterp xData yData] = interpCurvePDF(diffData,xq)
% Julian Rocha

% diffData is the raw diffusion coefficient/displacement data
% xq are query points

% Calculated the probability density function for the distribution
edges = xq;
binSpacing = edges(2)-edges(1);
histogram_counts = histcounts(diffData,(edges+binSpacing/2));
pdf_distribution = histogram_counts./(sum(histogram_counts));
pdf_distribution = pdf_distribution*(1/binSpacing);


xData = edges;
yData = zeros(1,length(xData));
yData(1,2:end) = pdf_distribution;

% Perform a B-spline interpolation of the data.
numPoints = length(diffData);
order = 25;
curveInterp = Bspline(order, numPoints, xData, yData);


% Query the curve at the query points in xq so that the output has a
% defined number of y values
xData = xq;
yData = curveInterp(xq);
end
