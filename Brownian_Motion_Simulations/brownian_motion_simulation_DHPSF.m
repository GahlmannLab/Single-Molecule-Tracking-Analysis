function x = brownian_motion_simulation_DHPSF ()
% brownian_motion_simulation_DHPSF simulates Brownian motion of diffusive
% molecules. This function produces images of the double-helix
% point-spread-function (DHPSF) based on the position of the molecule at a
% given time point. The images are generated to match experimentally
% acquired images. The Brownian motion is simulated to be confined to the
% volume of a cylinder (the approximate shape of a rod-shaped bacteria).

%Diffusion Coefficient to simulate, d
% for d = 1:1:15
for d = 10
    
    %The number of trajectories
    numTracks = 50;
    
    %l is half the total length of the cylinder (µm)
    l = 2.5;
    
    %r is the radius of the cylinder (µm)
    r = 0.4;
    
    %image size
    pixSize = 108; %nm
    imgSize = 2*(round(l*1000/pixSize)) + 51;
    
    
    %A look up table was created to sample the DHPSF at different z positions.
    %To save computational time, instead of explicitely simulating the DHPSF
    %for each point, the image with the closest z position in the look up table
    %is used.
    load('DHPSF_library.mat'); %simulated DHPSF library
    DHPSF_library = DHPSF_library;
    DHPSF_idx = DHPSF_idx;
    
    
    %calibration files from experimental data set (reflected channel)
    %same size as cell regions (97x97 pixels)
    load('DHPSF_darkImg.mat')
    DHPSF_darkImg = DHPSF_darkImg;
    
    load('DHPSF_readN.mat');
    DHPSF_readN = DHPSF_readN;
    
    load('DHPSF_gain.mat');
    DHPSF_gain = DHPSF_gain;
    
    
    %number of dimensions
    m = 3;
    %exposure time,s. This is set to match the experimental camera rate.
    exposure = 0.025;
    %dt is the short time step of the trajectory. This value must be much
    %lower than the exposure time.
    dt = 0.0000001;
    %track length = 6
    trackLength = 6;
    %add 2 to pad track
    n = (trackLength+2)*(exposure/dt);
    n = round(n);
    %total time of the simulation,s. Add 3 to the trackLength to add 3
    %frames in between tracks
    numFrames = numTracks*(trackLength+3);
    T = numFrames*exposure;
    
    %Mean and standard deviation for the photon distribution to sample
    %from.
    meanPhotons = 2000;
    photonsStd = 100;
    
    %separated time start, add a pad of 3 frames in between tracks
    time = (1:(trackLength+3):numTracks*(trackLength+3))';
    time = time*exposure;
    
    %tracksFinal is a cell array, with each cell containing 6 (trackLength)
    %x,y,z localizations for a trajectory
    tracksFinal = cell(numTracks,1);
    
    %Parfor loop can be used instead to decrease computational time
    parfor k = 1:numTracks
        %     for k = 1:numTracks
        
        %initialize temporary variables
        tempTrack = [];
        x = [];
        dx = [];
        xTemp = [];
        
        %Randomly set initial points of the track in the cylinder
        z0 = 2*r*rand(1,1) - r;
        y0(1,1) = 2*sqrt((r^2 - (abs(z0(1,1)))^2))*rand(1,1) - sqrt((r^2 - (abs(z0(1,1)))^2));
        x0 = 2*l*rand(1,1) - l;
        
        x = zeros(m,n);
        x(:,1) = [x0,y0,z0];
        
        %Compute the individual steps (short time-step of dt). Sample the
        %stepsize from a normal distribution.
        s = sqrt ( 2 * m * d * dt ) * randn ( 1, n - 1 );
        %  Choose a random direction for the molecule to move.
        if ( m == 1 )
            dx(1:m,1:n-1) = s(1:n-1);
        else
            a = randn ( m, n - 1 );
            v = s ./ sqrt ( sum ( a.^2 ) );
            b = spdiags ( v', 0, n-1, n-1 );
            dx(1:m,1:n-1) = a * b;
        end
        
        %Make sure the molecule remains in cylinder. If the molecule lands
        %outside of the cylinder, reflect it back inside the cylinder at a
        %random angle.
        for c = 1:n-1
            check = x(:,c) + dx(:,c);
            if abs(check(1,1)) > l
                if check(1,1) < 0
                    xTemp(1,1) = check(1,1)+l;
                    x(1,c+1) = -l - xTemp(1,1);
                else
                    xTemp(1,1) = check(1,1)-l;
                    x(1,c+1) = l - xTemp(1,1);
                end
            else
                x(1,c+1) = check(1,1);
            end
            if sqrt(check(2,1)^2 + check(3,1)^2) > r
                dTemp = sqrt(dx(2,c)^2 + dx(3,c)^2);
                yTemp = x(2,c);
                zTemp = x(3,c);
                pass = 0;
                while pass == 0
                    if sqrt(yTemp^2 + zTemp^2) > r
                        a = randn(2,1);
                        v = a*(dTemp/sqrt(sum(a.^2)));
                        yTemp = x(2,c) + v(1);
                        zTemp = x(3,c) + v(2);
                    else
                        pass = 1;
                    end
                end
                x(2,c+1) = yTemp;
                x(3,c+1) = zTemp;
            else
                x(2,c+1) = check(2,1);
                x(3,c+1) = check(3,1);
            end
        end
        
        tempTrack = x';
        
        %Inititalize variables
        photonsTot = [];
        xLoc = [];
        yLoc = [];
        zLoc = [];
        xLocRel = [];
        yLocRel = [];
        xMean = [];
        yMean = [];
        zMean = [];
        Frame = [];
        xLocRelTemp = [];
        yLocRelTemp = [];
        trackImg = cell(1,trackLength);
        
        
        %Sample from normal distribution for photon count, making sure
        %photon count is not less than 0
        photons = normrnd(meanPhotons, photonsStd,1);
        while photons < 0
            photons = normrnd(meanPhotons, photonsStd,1);
        end
        
        xCenter = imgSize/2;
        yCenter = imgSize/2;
        
        %idxTemp are the time indices of the trajectory (with short
        %time-steps) that correspond to the mid-point of the duration of
        %each of the frames
        idxTemp = single(1+(k-1)*((exposure/dt)/100):(exposure/dt):n);
        idxTemp = idxTemp(2:end-1);
        %tracks final are the localizations in the track located at the
        %mid-point of the duration of the frame
        tracksFinal{k} = tempTrack(idxTemp,:);
        
        if time(k,1) > T - ((n/(exposure/dt))-1)*exposure;
            tracksFinal{k} = [];
            tempTime = [];
        else
            tempTime = time(k,1):exposure:(time(k,1)+((n/(exposure/dt))-1)*exposure);
            tempTime = tempTime(2:end-1);
            
            for s = 1:length(idxTemp)
                %idxTemp2 are the indices of time points throughout the
                %duration of a frame. There are numAvg indices centered
                %around the indices in idxTemp. These are used to query
                %numAvg positions throughout the duration of the frame to
                %be averaged.
                numAvg = 50;
                idxTemp2 = single((idxTemp(s) - (exposure/dt)/2):(exposure/dt)/numAvg:(idxTemp(s) + (exposure/dt)/2));
                
                tracksTemp = tempTrack(idxTemp2,:);
                stepImg = zeros(imgSize);
                
                
                for t = 1:length(tracksTemp)
                    stepImgTemp = zeros(imgSize);
                    
                    xPix = ceil(xCenter + tracksTemp(t,1)*1000/pixSize);
                    yPix = ceil(yCenter + tracksTemp(t,2)*1000/pixSize);
                    
                    xRel = xCenter + tracksTemp(t,1)*1000/pixSize;
                    yRel = yCenter + tracksTemp(t,2)*1000/pixSize;
                    
                    x0 = (0.5 + (xRel - xPix));
                    y0 = (0.5 + (yRel - yPix));
                    
                    %Based on the position of the molecule, obtain an image
                    %of the DHPSF from the DHPSF_library.
                    temp = abs(DHPSF_idx - tracksTemp(t,3));
                    [temp2 tempIdx] = min(temp);
                    tempImg = DHPSF_library{tempIdx};
                    tempImg = fftshift(fft2(tempImg));
                    [xF,yF] = meshgrid(-15:15,-15:15);
                    tempImg = tempImg.*exp(-1i*2*pi.*(xF*x0+yF*y0)/31);
                    tempImg = ifft2(tempImg);
                    tempImg = abs(tempImg);
                    
                    xLocRelTemp = (0.5 + (xRel - xPix))*pixSize;
                    yLocRelTemp = (0.5 + (yRel - yPix))*pixSize;
                    
                    %stepImgTemp is the temporary image representing the
                    %DHPSF at that specific time point. stepImg is the sum
                    %of the DHPSF images at numAvg time points
                    stepImgTemp(yPix-15:yPix+15,xPix-15:xPix+15) = tempImg;
                    stepImg = stepImg + stepImgTemp;
                end
                
                %Normalize the image so that the sum of all pixels is equal
                %to the value in photons
                stepImg = stepImg/(sum(sum(stepImg)));
                stepImg = stepImg*photons;
                
                %trackImg contains the images for the 6 (trackLength)
                %frames in the trajectory. At this point noise, camera
                %background, etc, has not been added to the image
                trackImg{s} = stepImg;
                
                %The mean localizations are the mean position of the
                %positions at the times indexed by idxTemp2
                xMean = vertcat(xMean,mean(tracksTemp(:,1)));
                yMean = vertcat(yMean,mean(tracksTemp(:,2)));
                zMean = vertcat(zMean,mean(tracksTemp(:,3)));
                
                
                xLocRel = vertcat(xLocRel,xLocRelTemp);
                yLocRel = vertcat(yLocRel,yLocRelTemp);
            end
            
            %all x,y,z localizations and photon numbers are stored in these
            %variables
            photonsTot = vertcat(photonsTot,photons);
            xLoc = vertcat(xLoc,tracksFinal{k}(:,1));
            yLoc = vertcat(yLoc,tracksFinal{k}(:,2));
            zLoc = vertcat(zLoc,tracksFinal{k}(:,3));
            
        end
        
        
        %fiberData is a structure array containing all important
        %information for the trajectories
        Frame = vertcat(Frame,tempTime');
        fiberData(k).FrameRange = [1 T/exposure];
        fiberData(k).xLoc = xLoc*10^3;
        fiberData(k).yLoc = yLoc*10^3;
        fiberData(k).zLoc = zLoc*10^3;
        fiberData(k).xLocRel = xLocRel;
        fiberData(k).yLocRel = yLocRel;
        fiberData(k).zLocRel = zLoc*10^3;
        fiberData(k).xMean = xMean*10^3;
        fiberData(k).yMean = yMean*10^3;
        fiberData(k).zMean = zMean*10^3;
        fiberData(k).photons = photonsTot;
        fiberData(k).Frame = round((Frame/exposure));
        fiberData(k).trackImg = trackImg;
        
        k
    end
    clear s v errorX errorY errorZ DHPSF_library DHPSF_idx a b x xTemp stepImg stepImgTemp
    
    totalImg = cell(numFrames,1);
    totalBkgndImg = cell(numFrames,1);
    
    trackImgTot = [];
    frameTot = [];
    for b = 1:numTracks
        trackImgTot = vertcat(trackImgTot,fiberData(b).trackImg');
        frameTot= vertcat(frameTot,fiberData(b).Frame);
    end
    
    for c = 1:numFrames
        [goodStep goodIdx] = find(frameTot == c);
        temp1 = zeros(imgSize);
        if ~isempty(goodStep)
            for e = 1:length(goodStep)
                if isempty(trackImgTot{goodStep(e)})
                    trackImgTot{goodStep(e)} = zeros(imgSize);
                    
                end
                temp1 = temp1 + trackImgTot{goodStep(e)};   
            end
        end
        regionTemp{c} = temp1;
    end
    
    %divide the DHPSF_darkImg by the gain of the camera to get an image
    %in units of photons (it is loaded in as an image in units of
    %counts).
    DHPSF_darkImg_photons = DHPSF_darkImg./DHPSF_gain;
    
    
    for b = 1:numFrames
        %add poisson noise to background, 13.5 photons determined
        %experimentally
        DHPSF_bkgndImg = poissrnd(13.5,imgSize,imgSize);
        
        %add poisson noise to DHPSF image
        temp = regionTemp{b};
        temp = poissrnd(temp);
        
        %calculate camera read noise
        readNoise = normrnd(zeros(imgSize),DHPSF_readN);
        readNoise = readNoise./DHPSF_gain;
        
        %add background with poisson noise, DHPSF image with poisson noise,
        %camera read noise, and the darkImg
        temp = temp + DHPSF_bkgndImg + readNoise + DHPSF_darkImg_photons;
        temp2 = DHPSF_bkgndImg + readNoise + DHPSF_darkImg_photons;
        
        %multiply by camera gain to get image in units of camera counts
        temp = temp.*DHPSF_gain;
        temp2 = temp2.*DHPSF_gain;
        
        temp = uint16(temp);
        temp2 = uint16(temp2);
        
        %totalImg contains all data frames
        totalImg{b} = temp;
        %totalImg contains the simulated background for all frames
        totalBkgndImg{b} = temp2;
        
    end
    
    
    clear regionTemp temp temp2 trackImgTot
    
    DC = cell(200,1);
    %create 200 frames of dark images
    for a = 1:200
        readNoise = normrnd(zeros(imgSize),DHPSF_readN);
        readNoise = readNoise./DHPSF_gain;
        %multiply by camera gain to get an image in units of counts
        DC{a} = uint16((DHPSF_darkImg_photons + readNoise).*DHPSF_gain);
    end
    
    %name the save file
    filename = ['d', num2str(d),'_r', num2str(r*1000),'nm_dt',...
        num2str(exposure*1000),'ms_',num2str(numTracks),'tracks_DHPSF2'];
    
    
    %remove images from fiberData so that the file can be saved properly
    fiberData = rmfield(fiberData,'trackImg');
    %save file
    
    filename = strrep(filename, '.', '-');
    save([filename '.mat'],'-v7.3');
    
end