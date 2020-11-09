function x = brownian_motion_simulation_DHPSF_membrane_main ()
% Julian Rocha
% brownian_motion_simulation_DHPSF simulates Brownian motion of diffusive
% molecules. This function produces images of the double-helix
% point-spread-function (DHPSF) based on the position of the molecule at a
% given time point. The images are generated to match experimentally
% acquired images. The Brownian motion is simulated to be confined to the
% volume of a cylinder (the approximate shape of a rod-shaped bacteria).
for iter = 1:5
%Diffusion Coefficient to simulate, d
% for d = 1:1:15
% for d = 0.5
%for cytoF = 0:0.2:1
    for cytoF = 0
        dM = 0.2;
        dC = 0.2;
    %The total number of trajectories
    numTracks = 1000;
    
    %%Additions and modifications by Alecia Achimovich 5/13/2020
    %Cytosolic fraction of tracks
    %cytoF = 0;
        
    %%------End modifications by Alecia Achimovich -----
    
    %l is half the total length of the cylinder (µm)
    %Modified length of cylinder so that the length of the cell including
    %the spherical ends is 5 um. Length of cylinder prior to modification =
    %2.5 um. Alecia Achimovich 5/20/2020
    l = 2.1;
    
    %r is the radius of the cylinder (µm)
    r = 0.4;
    
    %image size
    %Modified image size to include hemispherical ends. Previous equation:
    %imgSize = 2*(round(l*1000/pixSize))+51;
    pixSize = 108; %nm
    imgSize = 2*(round((l+r)*1000/pixSize)) + 51;
    
    
    %A look up table was created to sample the DHPSF at different z positions.
    %To save computational time, instead of explicitely simulating the DHPSF
    %for each point, the image with the closest z position in the look up table
    %is used.
    load('DHPSF_library.mat'); %simulated DHPSF library
    DHPSF_library = DHPSF_library;
    DHPSF_idx = DHPSF_idx;
    
    
    %calibration files from experimental data set (reflected channel)
    %same size as cell regions (97x97 pixels)
    load('DHPSF_darkImg_2.mat')
    DHPSF_darkImg = DHPSF_darkImg;
    
    load('DHPSF_readN_2.mat');
    DHPSF_readN = DHPSF_readN;
    
    load('DHPSF_gain_2.mat');
    DHPSF_gain = DHPSF_gain;
    
    
    %number of dimensions
    m = 3;
    %exposure time,s. This is set to match the experimental camera rate.
    exposure = 0.025;
    %dt is the short time step of the trajectory. This value must be much
    %lower than the exposure time.
    dt = 0.0000001;
    %Chen algorithm integration time
    %track length = 6
    %Modification by Alecia Achimovich - to more accurately reflect JF549
    %fluorophore
    trackLength = 10;
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
    
    tracksFine = cell(numTracks,1);
    kc = 1;
    ks = 1;
    
    %Parfor loop can be used instead to decrease computational time
    for k = 1:numTracks
        %     for k = 1:numTracks
        
        %initialize temporary variables
        tempTrack = [];
        x = [];
        dx = [];
        xTemp = [];
        %Randomly select whether fluorophore is membrane bound or in
        %cytosol
        corm = rand(1,1);       
        if corm <= cytoF
        %Randomly set initial points of the track in the cylinder
        %z0 = 2*r*rand(1,1) - r;
        %y0(1,1) = 2*sqrt((r^2 - (abs(z0(1,1)))^2))*rand(1,1) - sqrt((r^2 - (abs(z0(1,1)))^2));
        % --------------Below modified by Ting Yan------------------
        rnd_r = r*(rand(1,1)^(1/2));
        rnd_theta = 2*pi*rand(1,1);
        z0 = rnd_r*sin(rnd_theta);
        y0 = rnd_r*cos(rnd_theta);
        % --------------modification done------------------
        
        x0 = 2*l*rand(1,1) - l;
        
        x = zeros(m,n);
        x(:,1) = [x0,y0,z0];
        
        %Compute the individual steps (short time-step of dt). Sample the
        %stepsize from a normal distribution.
        s1 = sqrt ( 2 * m * dC * dt ) * randn ( 1, n - 1 );
        %  Choose a random direction for the molecule to move.
        if ( m == 1 )
            dx(1:m,1:n-1) = s1(1:n-1);
        else
            a = randn ( m, n - 1 );
            v = s1 ./ sqrt ( sum ( a.^2 ) );
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
                 %--------------below modified by Ting Yan----------
%                 syms ycrx zcrx real
%                 [Sy, Sz] = solve(ycrx^2 + zcrx^2 == r^2, (ycrx - x(2,c))/dx(2,c) ==(zcrx - x(3,c))/dx(3,c), [ycrx, zcrx], 'Real', true);
%                 y_cross = double(Sy(sign(Sy-x(2,c)) == sign(dx(2,c)))); % the point where the molecule hits the boundary
%                 z_cross = double(Sz(sign(Sy-x(2,c)) == sign(dx(2,c))));
%                 if isempty(y_cross)
%                     c
%                     x(:,c)
%                     dx(:,c)
%                     error
%                 end
                
%                 reflected_d = dTemp - sqrt((y_cross - x(2,c))^2 +(z_cross-x(3,c))^2 ); % the length that is reflected into the cylinder
%                 a = randn(2,1);
%                 v = a*(reflected_d/sqrt(sum(a.^2)));
%                 yTemp = y_cross + v(1); % new position y
%                 zTemp = z_cross + v(2); % new position z
                
%                 while pass == 0
%                     if sqrt(yTemp^2 + zTemp^2)> r
%                         a = randn(2,1);
%                         v = a*(reflected_d/sqrt(sum(a.^2)));
%                         yTemp = y_cross + v(1);
%                         zTemp = z_cross + v(2);
%                     else
%                         pass = 1;
%                     end
%                 end
                %--------------modification done----------
              
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
        tracksFine{k} = tempTrack;
        tracksFineCyl{kc} = tempTrack;
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
        %idxTemp = single(1+(k-1)*((exposure/dt)/100):(exposure/dt):n);
        idxTemp = single(round((exposure/dt)/2):exposure/dt:n); %Modification by Ting
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
        kc = kc+1;
        elseif corm > cytoF
%%Section below written and modified by Alecia Achimovich. For simulation of diffusion on the surface of the cell.     

% Set-up probability of finding component on cylinder of sphere of cell

%l is half full length of cell. h = full
h = 2*l;
SA_sphere = 4*pi*r^2;
SA_cylinder = 2*pi*r*h;
p_cylinder = SA_cylinder/(SA_cylinder+SA_sphere);
p_sphere = SA_sphere/(SA_cylinder+SA_sphere);

deltaT = [];
deltaST = [];
deltayzT = [];
mn=3;
   
        %     for k = 1:numTracks
        %initialize temporary variables
        tempTrack = [];
        tempTrackSphere = [];
        x = [];
        dx1 = [];
        xTemp = [];
        xcap = [];
        
        %Randomly determine if  point is on cylinder or end caps of cell. 
        CorS = rand(1,1);
               
        %Randomly set initial points of the track in the cylinder
        z0 = 2*r*rand(1,1) - r;
       
        %If initial point is in cylinder, y is confined to outer boundary. 
        if CorS > p_sphere
        x0 = 2*l*rand(1,1) - l;
        %y0(1,1) = sqrt((r^2 - (abs(z0(1,1)))^2));
        initTheta = 2*pi*rand(1,1);
        z0 = r*sin(initTheta);
        y0 = r*cos(initTheta);
        %Randomly determine if y is positive or negative
%         Coin = binornd(1,0.5);         
%             if Coin == 0
%                 y0 = y0;
%             elseif Coin == 1
%                 y0 = y0*(-1);
%             end
        %If initial point is in end cap, x is determined by y & z.     
        elseif CorS <= p_sphere
         %y0(1,1) = 2*sqrt((r^2 - (abs(z0(1,1)))^2))*rand(1,1) - sqrt((r^2 - (abs(z0(1,1)))^2));
         %xcyl(1,1) = sqrt((r^2 - (abs(z0(1,1)))^2)-(abs(y0(1,1))^2));
        s_angle = 2*pi*rand(1,1);
        %azimuth_angle = pi*rand(1,1);
        z0 = 2*r*rand(1,1) - r;
        xcyl = sqrt(r^2-z0^2)*cos(s_angle);
        y0 = sqrt(r^2-z0^2)*sin(s_angle);
%         xcyl(1,1) = r*(sin(elevation_angle))*cos(azimuth_angle);
%         y0 = r*sin(elevation_angle)*(sin(azimuth_angle));
%         z0 = r*cos(elevation_angle);
        %Randomly determine if x is positive or negative
%         Coin = binornd(1,0.5);         
%             if Coin == 0
%                 xcyl = xcyl;
%             elseif Coin == 1
%                 xcyl = xcyl*(-1);
%             end
            xcap(1,1) = xcyl;
            %Add end cap localization to end of cells.
            if xcyl < 0
                x0 = (-1*l)+xcyl;
            elseif xcyl > 0
                x0 = l+xcyl;
            end
        end
            
        x = zeros(m,n);
        x(:,1) = [x0,y0,z0];
        allLocs(:,k) = [x0,y0,z0];
        %Compute the individual steps (short time-step of dt). Sample the
        %stepsize from a normal distribution.
        s1 = sqrt ( 2 * mn * dM * dt ) * abs(randn ( 1, n - 1 )); 
       %1D stepsize for 2D diffusion along cylinder
        s1(2,:) = sqrt( 2 * dM * dt ) * randn ( 1, n - 1 );
        s1(3,:) = sqrt ( 2 * dM * dt ) * randn ( 1, n - 1 );
        %  Choose a random direction for the molecule to move.
        
        %Create an azimuth angle for each point in the sphere. 
        %These values are uniformly distributed
        %Points that are distributed along length of cylinder. Curvature is
        %along y-axis.
        z2 = [];
        z2(1,1) = x(3,1);
        x2 = [];
        x2(1,1) = x(1,1);
        y2 = [];
        y2(1,1) = x(2,1);
        for ll = 1:n-1
            
            if abs(x2(1,ll)) <= l % on the surface of cylinder                 
                [x2(1,ll+1),y2(1,ll+1),z2(1,ll+1)] = CylMemDiff(r,l,s1(2,ll),s1(3,ll),x2(1,ll),y2(1,ll),z2(1,ll));                
                
                %Make sure that molecule doesn't diffuse away from membrane.                
              
                %Points on endcaps.
            xcap = [];
            elseif abs(x2(1,ll)) > l
                phi = 2*pi*rand(1,1);
                
                [x2(1,ll+1),xcap(1,ll+1),y2(1,ll+1),z2(1,ll+1)] = EndCapMemDiff(r,l,phi,s1(1,ll),z2(1,ll),x2(1,ll),y2(1,ll));
%                 if x2(1,ll) < 0 && x2(1,ll+1) > 0 || x2(1,ll) > 0 && x2(1,ll+1) < 0
%                     
%                     dxover = abs(x2(1,ll+1))-l;
%                         if x2(1,ll) > 0 && x2(1,ll+1) < 0                            
%                             x2(1,ll+1) = l-dxover;
%                         elseif x2(1,ll) < 0 && x2(1,ll+1) > 0
%                             x2(1,ll+1) = -l+dxover;
%                         end
%                 end
            end
        end   
          
        
    %end
       
        tempTrack(1,:) = x2;
        tempTrack(2,:) = y2;
        tempTrack(3,:) = z2;
        tempTrack = tempTrack';
        tracksFine{k} = tempTrack;
        tracksFineSphere{ks} = tempTrack;
%         delta = [];
%         deltaS = [];
%         deltayz = [];
%         rs = 1;
%         rt = 1;
%         for rl = 1:length(tempTrack)-1
%             if abs(tempTrack(rl,1)) < l
%             delta(rs,1) = tempTrack(rs+1,1)-tempTrack(rs,1);
%             delta(rs,2) = tempTrack(rs+1,2)-tempTrack(rs,2);
%             delta(rs,3) = tempTrack(rs+1,3)-tempTrack(rs,3);
%             deltayz(rs,1) = delta(rs,1);
%             deltayz(rs,2) = sqrt(delta(rs,2)^2+delta(rs,3)^2);
%                 if delta(rs,3) < 0
%                     deltayz(rs,2) = -deltayz(rs,2);
%                 end
%             rs = rs + 1;
%             elseif abs(tempTrack(rl,1)) > l
%             deltaS(rt,1) = tempTrack(rt+1,1)-tempTrack(rt,1);
%             deltaS(rt,2) = tempTrack(rt+1,2)-tempTrack(rt,2);
%             deltaS(rt,3) = tempTrack(rt+1,3)-tempTrack(rt,3);
%             rt = rt+1;
%             end
%         end
%       deltaT = [deltaT;delta];
%       deltaST = [deltaST;deltaS];
%       deltayzT = [deltayzT;deltayz];
%         
        %tempTrackSphere(1,:) = xcap;
%         tempTrackSphere(2,:) = y2;
%         tempTrackSphere(3,:) = z2;
%         tempTrackSphere = tempTrackSphere';
%         tracksFineSphere{k} = tempTrackSphere;
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
        %idxTemp = single(1+(k-1)*((exposure/dt)/100):(exposure/dt):n);
        idxTemp = single(round((exposure/dt)/2):exposure/dt:n); %Modification by Ting Yan
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
        ks = ks+1;
        end
    end
%     for jk = 1:length(tracksFine)
%     for ll = 1:length(tracksFine{jk})-1
%     delta(ll,1) = tracksFine{jk}(ll+1,1)-tracksFine{jk}(ll,1);
%     delta(ll,2) = tracksFine{jk}(ll+1,2)-tracksFine{jk}(ll,2);
%     delta(ll,3) = tracksFine{jk}(ll+1,3)-tracksFine{jk}(ll,3);   
%     end
%     deltaT = [deltaT;delta];
%     delta = [];
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
    filename = ['cytoF', num2str(cytoF),'_dM',num2str(dM),'_dC',num2str(dC),'_r', num2str(r*1000),'nm_dt',...
        num2str(exposure*1000),'ms_',num2str(numTracks),'tracks_DHPSF2_',num2str(iter)];    
    
   %remove images from fiberData so that the file can be saved properly
    fiberData = rmfield(fiberData,'trackImg');
    %save file
    tic
    filename = strrep(filename, '.', '-');
    save([filename '.mat'],'-v7.3');
    toc
    clear fiberData
    end
end
end
