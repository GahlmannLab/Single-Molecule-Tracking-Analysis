% Copyright (c)2013, The Board of Trustees of The Leland Stanford Junior
% University. All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
% Neither the name of the Leland Stanford Junior University nor the names 
% of its contributors may be used to endorse or promote products derived 
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [laser_x_nm, laser_y_nm ,sigma_x_nm, sigma_y_nm, theta, peakIntensity, waist]...
    = f_findGaussianLaserProfile...
    (bkgndImg_avg, FOVmask, nmPerPixel, powerAtObjective, ROI)
%UNTITLED2 Summary of this function goes here
%   Estimate Gaussian laser background 


[xIdx yIdx] = meshgrid(1:size(bkgndImg_avg,2), 1:size(bkgndImg_avg,1));

bkgndImg_avg_masked =  bkgndImg_avg.*FOVmask;

bkgndFit = [max( bkgndImg_avg_masked( bkgndImg_avg_masked>0))- min( bkgndImg_avg_masked( bkgndImg_avg_masked>0)), ...
    size( bkgndImg_avg,2)/2, size( bkgndImg_avg,1)/2, ...
    size( bkgndImg_avg,2)/4, size( bkgndImg_avg,1)/4, 0, ...
    min( bkgndImg_avg_masked( bkgndImg_avg_masked>0))];

% Fit with lsqnonlin
lowerBound = [0, 1, 1, 1, 1, -360, 0];
upperBound = [2*(max( bkgndImg_avg_masked( bkgndImg_avg_masked>0))- min( bkgndImg_avg_masked( bkgndImg_avg_masked>0))),...
    size( bkgndImg_avg,2), size( bkgndImg_avg,1),...
    2*size( bkgndImg_avg,2), 2*size( bkgndImg_avg,1),...
    360,...
    2*min( bkgndImg_avg_masked( bkgndImg_avg_masked>0))];

bkgndFit = lsqnonlin(@(x) singleGaussianRotOffset(x, bkgndImg_avg,xIdx,yIdx,FOVmask),...
    bkgndFit,lowerBound,upperBound);

laser_x = bkgndFit(2);
laser_y = bkgndFit(3);
sigma_x = bkgndFit(4);
sigma_y = bkgndFit(5);
theta = bkgndFit(6)*pi/180;
A = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
B = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
C = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;

bkgndFit_Image = bkgndFit(1)*exp( -(A*(xIdx-bkgndFit(2)).^2 +...
    2*B*(xIdx-bkgndFit(2)).*(yIdx-bkgndFit(3)) +...
    C*(yIdx-bkgndFit(3)).^2)) +bkgndFit(7);

hFitProfileFig = figure;
subplot(1,2,1)
imagesc( bkgndImg_avg.*FOVmask);
axis square, colorbar
title('Average Background Image (masked)');
subplot(1,2,2)
imagesc(bkgndFit_Image.*FOVmask);
axis square, colorbar
title('Fit to Average Background');

% shift coordinates relative to entire dataset (not just ROI)
laser_x_nm =  (laser_x + ROI(1)-1) * nmPerPixel;            % laser center position
laser_y_nm =  (laser_y + ROI(2)-1) * nmPerPixel;
sigma_x_nm = sigma_x * nmPerPixel;    % Gaussian Beam waist (1/e^2) is 2*sigma
sigma_y_nm = sigma_y * nmPerPixel;

waist = 2*mean([sigma_x, sigma_y]) * nmPerPixel / 1000; % in units of microns

% Assuming a radially symmetric Gaussian intensity distribution
peakIntensity = 4*(powerAtObjective)/(2*pi*waist^2)...
    * (10^4)^2;  % in units of Watts/cm^2


end

function err = singleGaussianRotOffset(par,Zdata,ii,jj,mask)
sigma_x = par(4);
sigma_y = par(5);
theta = par(6)*pi/180;

A = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
B = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2;
C = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;

Zmodel = par(1)*exp( -(A*(ii-par(2)).^2 + 2*B*(ii-par(2)).*(jj-par(3)) + C*(jj-par(3)).^2)) ...
    +par(7);
%Zmodel = par(1).*exp( -((ii-par(2)).^2+(jj-par(3)).^2) / (2*par(4).^2)) ...
%    +bkgnd;

% lower = Zdata>threshold(1);
% upper = Zdata<threshold(2);
% good = lower & upper;

err = reshape(Zmodel(mask)-Zdata(mask),1,[]);
end
