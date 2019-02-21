function [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dataFile,medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt)
% currIdx = find(frames==frames(c)); % current index within 'frames'
%after this line, the darkcounts/ offset is taken out
if isDcimg
    [framedata,~]= dcimgmatlab(frames(c)-1, dataFile);
    framedatatrans = transpose (framedata);
%     totalframes = double(totalframes);
    data = double(framedatatrans)-darkAvg;
else
    data = double(imread(dataFile,frames(c))) - darkAvg;
end

%             dataLaser = data;
%after this line, data considered is strictly data in the square
%ROI
data = data(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
%After this line, data is officially in photons not counts
%data = double(data*conversionFactor);

if ~isequal(size(data),size(gain))
    gain = gain(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
end
data = double(data./gain);

if nhaData && usePolyROI
    data=data.*(FOVmask1);
end

if usePolyROI
    %data ring is the measure in photons of data in an external
    %ring of the polygon ROI that was drawn (the ring is about 40
    %pixels thick) the median is then applied as the "data to all
    %of the pixels outside the polygon ROI
    dataRing=data(FOVmask&~FOVmask1);
    data(~FOVmask)=median(dataRing);
    % the below can smooth the data near the edge of the FOV -
    % would be useful if using wavelet bg subtraction, e.g.
    %             for i = 1:10
    %             dataBlur=imfilter(data,gaussFilt,'replicate');
    %             data(~FOVmask)=dataBlur(~FOVmask);
    %             end
    %             clear dataBlur
end

if medianFilter   
    if mod(frames(c),interpVal)==0
        bkgndImg = bkgndImgvals{idivide(int32(frames(c)), int32(interpVal)), 1};
    else
        bkgndImg = bkgndImgvals{idivide(int32(frames(c)), int32(interpVal))+1, 1};
    end
    if medianBlurSigma~=0
        bkgndImg(~FOVmask)=median(reshape(bkgndImg(~FOVmask1&FOVmask),[],1));
        bkgndImg = imfilter(bkgndImg,medBlurFilt,'replicate');
    end
    bkgndImg = bkgndImg./gain;    
elseif nhaData
    % use median of all pixels as BG estimate
    bkgndImg = median(data(:)).*ones(size(data));
else
    bkgndImg = f_waveletBackground(data);
    %                 bkgndImgLaser = f_waveletBackground(dataLaser);
end

if usePolyROI && medianFilter % otherwise the displayed image is wrong
%     gain1 = gain(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
    %data is multiplied by gain so bkgndImg will again be in counts
    %not photons yet in and out of the ROI
%     bkgndImg(~FOVmask) = data(~FOVmask)./gain1(~FOVmask);
    bkgndImg(~FOVmask) = data(~FOVmask)./gain(~FOVmask);

end
%             bkgndImgTotalLaser = bkgndImgTotalLaser + bkgndImgLaser;
bkgndImg = double(bkgndImg);
data = data - bkgndImg;
end

