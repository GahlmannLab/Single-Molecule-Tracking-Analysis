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

% f_medianFilter can be called within f_calSMidentification or f_fitSMs. It
% takes a set of frames and generates a rolling array named dataWindow
% where the first dimension is the stack of frames within the window, and
% the second and third dimensions are the Y and X of the image.
%
% The median of the dataWindow array along the first dimension then be used
% to generate a median image of the data.
%
% In practice, this image is generally noisy and should be smoothed. (This
% may not be the case if using windows >> 100 frames in width.) This is
% done in the parent functions. A Gaussian kernel with 15-pixel sigma seems
% to work so far. Could use a smaller kernel, but this one seems to be
% safer: noise is more mitigated, and faint images of the DHPSF that appear
% when using 100-frame windows (if not larger windows) are smoothed out.
%
% It is possible to weight the frames by their mean brightness, but frame
% to frame brightness variability was low in the test data analyzed and so
% this option is off by default.
%
% Note that 'frames' array likely does not start at 1, and is unlikely to
% be continuous, especially for interleaved multicolor data. Thus, we
% generate a subindex within the 'frames' array called 'winIdx.' winIdx is
% a list of the frames within 'frames' that are to be used in the window.
%
% Loading images from a .tif is a slow process, so dataWindow is changed,
% not recreated each time. To add and remove frames from dataWindow,
% dataWindow is circularly permuted and only new frames (near the end, if
% moving forward, or near the start, if moving backward, as for
% thresholding) are loaded. This function currently supports going backward
% or forward, with arbitrary spacing. This is done by remembering the last
% frame in the parent function, and then using the spacing between the last
% and current frames to update the background image.
%
% Interpolation is done within f_fitSMs and f_calSMidentification by only
% calling this function every N frames, and using the last-generated
% background image in the meantime. (The 'N frames' is defined within the
% 'frames' array: if we specify 5-frame interpolation and the frame list is
% [1, 2, 11, 12, 21, 22, 31, 32, 41, 42, 51,...], a median background image
% would generated on frame 1 and used for frames 1, 2, 11, and 12; then
% another median background image would be generated centered on 21 and
% used until frame 42, etc.) This could be switched to 'true'
% interpolation, interpolating pixel-by-pixel within the background images,
% but this method is used to keep memory costs low (no need to store more
% than one background image)

function [bkgndImg,dataWindow] = f_medianFilter_Sim(totalImg, darkAvg, ROI, frames, currFrame, windowSize, dataWindow,lastFrame)
    

    % do not use a more complicated iterating process or weight frames:    
    simpleMedian = 1;
    weightMedian = 0;
    numIter = 1;
    %
    % pass in dataWindow as nan first time, then keep the dataWindow array in
    % the parent function to feed back in.

    % windowSize is total size front and back. Even if specified as even,
    % will end up being odd (include central frame in addition to front and
    % back)
    
    % index of current frame within 'frames' list
    currIdx = find(frames==currFrame); 
    lastIdx = find(frames==lastFrame);
    
    % index of all frames within the window within the 'frames' list
    winIdx = currIdx-floor(windowSize/2):currIdx+floor(windowSize/2);
    
    % keep nans as placeholders when window goes beyond edges of 'frames'
    winIdx(winIdx<1|winIdx>length(frames))=nan;
    
    % generate dataWindow x imageHeight x imageWidth array for median
    
    
    if isnan(lastFrame)
        %% no usable previous frame, we must generate window de novo
        for bgFr = winIdx
            relFr = find(winIdx==bgFr); % index within window
            if isnan(bgFr)
                dataWindow(relFr,:,:) = nan;
            else 
                dataFr = frames(bgFr);
                newData = double(totalImg{dataFr})-darkAvg;
                newData = newData(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
                dataWindow(relFr,:,:) = newData;
            end
        end
    else
        %% add only new frames at the top, remove the ones at the bottom
        stepSize = currIdx-lastIdx;
        dataWindow = circshift(dataWindow,-stepSize);
         % the dataset moves in the opposite direction of the frame.
         % analogy/mnemonic: if the car moves forward, the road moves backward
         %
         % newWin is the set of frames within winIdx (that is, the set of
         % frames within the window) that are new and will have to be
         % replaced in dataWindow after permuting the elements of
         % dataWindow.
         %
         % this supports abs(stepSize) > 1 
         %
        if stepSize > 0 % frames moving forward (fitting)
            newWin = length(winIdx):-1:length(winIdx)-(stepSize-1);
        elseif stepSize < 0 % frames moving backward (threshold)
            newWin = 1:-stepSize;
        end
            for j = newWin
                if isnan(winIdx(j))
                    dataWindow(j,:,:) = nan;
                else
                    dataFr = frames(winIdx(j));
                    newData = double(totalImg{dataFr})-darkAvg;
                    newData = newData(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
                    dataWindow(j,:,:) = newData;
                end
            end
    end
    % note that median values are not integer due to subtraction of
    % darkAvg; however, the possible values in each pixel should be integer
    % + a double that is constant for the entire dataset
    if simpleMedian
        %% default, likely to stay this way
        bkgndImg = squeeze(nanmedian(dataWindow,1));
    elseif weightMedian 
        %% alternative, untested method
        % weight as in Hoogendoorn et al., Sci Rep 2014
        % seems to not be necessary: true background values are ~constant,
        % while the mean value jitters based on #SMs in frame (not from
        % changes that would alter the true bg, like laser power
        % fluctuation, etc.) Thus, scaling based on mean value without
        % somehow removing the SM constribution to the mean is
        % inappropriate for our data.
        meanVal=mean(mean(dataWindow(1,:,:),3),2);
        meanFrame = mean(reshape(dataWindow(1+floor(windowSize/2),:,:),[],1));
        bkgndImg=squeeze(nanmedian(dataWindow.*repmat(meanVal,[1,size(dataWindow,2),size(dataWindow,3)]),1));
        bkgndImg = bkgndImg.*meanFrame;
    elseif ~simpleMedian && numIter>0
        %% iteratively reduce background in hot pixels: not fully tested
        % seems to be limited, can't always fully remove image of DHPSF.
        % also can be slow (> 1 second / iteration)
%         medval = nanmedian(dataWindow(:,189,175));
        zeroT = tic;
        rangeFactor = 1; % 2.355
%         figure; hist(dataWindow(:,189,175),0:5:120); title(['0th iteration distribution: 80+- 50 frames, xy = 175,189, med=' num2str(medval) 'cts']);
%         hold on; plot([medval medval],[0,30],'r');
        for i = 1:numIter
        dwStd = nanstd(dataWindow,0,1); % 1 x height x width
        dwMed = nanmedian(dataWindow,1);
        dwRange = cat(1,dwMed - dwStd*rangeFactor, dwMed + dwStd*rangeFactor);
        dwOK = dataWindow >= repmat(dwRange(1,:,:),[length(winIdx),1,1]) &...
               dataWindow <= repmat(dwRange(2,:,:),[length(winIdx),1,1]);
        dataWindow(~dwOK) = nan;
        %dwCorrMed = squeeze(nanmedian(dataWindow,1));
        toc(zeroT)
        end
        
        bkgndImg = squeeze(nanmedian(dataWindow,1));
    end
    
end