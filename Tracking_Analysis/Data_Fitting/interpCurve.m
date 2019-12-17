function [curveInterp xData yData] = interpCurve(diffData,xq)
% Julian Rocha

% diffData is the raw diffusion coefficient/displacement data
% xq are the query points

curveInterp = [];
curveTemp = [];
curveTempX = [];

% Calculate eCDF (Empirical cumulative distribution funciton) of the 
% distribtuion
[curveTemp, curveTempX, flo, fup] =  ecdf(diffData,'alpha',0.32,'bounds','on');
curveTempXi = curveTempX;

fIdx = (1:length(flo))';

% Pad curvesTemp on both sides so that the Bspline interpolant can be
% queried outside of the range of the ecdf curve
curveTempXi = curveTempX;
[curveTempX curveTemp] = padArrayInterp(curveTempXi, curveTemp, xq);

% Interpolate ecdf so that all curves have the same number of points,
% necessary for scatteredInterpolant function later
order = 10;
numPoints = 10000;
curveInterp = Bspline(order,numPoints,curveTempX,curveTemp);


% Query the curve at the query points in xq so that the output has a
% defined number of y values
xData = xq;
yData = curveInterp(xq);
end
