nmPerPixel = 108;
channel = 'g';

[dataFile, dataPath] = uigetfile({'*.dcimg;*.tif;*.tiff','Standard image types';'*.*','All Files'},'MultiSelect','on','Open image stack for data processing');

[darkFile, darkPath] = uigetfile({'*.dcimg;*.tif;*.tiff','Standard image types';'*.*','All Files'},'Open image stack with dark counts (same parameters as calibration)',dataPath);
   

temp = inputdlg({'What was the laser power at the objective? (in mW)'},...
                'Input laser power',...
                1,...
                {'0'});
   powerAtObjective = str2double(temp{1})/1000;
   
   
   lastDir=dataPath;
        dcimgfile = fullfile(dataPath, dataFile);
        dcimgfile = strrep(dcimgfile, '\', '\\');
        [framedata,totalframes]= dcimgmatlab(0, dcimgfile);
        totalframes = double(totalframes);
        [imgHeight, imgWidth] = size(transpose(framedata));
        numFrames = double(totalframes);
        numFramesTotal = numFrames;
        numFiles =1;
        
        dcimgfile_D = fullfile(darkPath, darkFile);
            dcimgfile_D = strrep(dcimgfile_D, '\', '\\');
            [framedata_dark,numDarkFrames]= dcimgmatlab(0, dcimgfile_D);
            [darkHeight, darkWidth] = size(transpose(framedata_dark));
            darkAvg = zeros(darkHeight, darkWidth);
            for frame = 0:numDarkFrames-1
                [framedata_dark,numDarkFrames]= dcimgmatlab(frame, dcimgfile_D);
                framedatatrans = transpose (framedata_dark);
                darkAvg = darkAvg + double(framedatatrans);
            end
            darkAvgLaser = darkAvg./double(numDarkFrames);
   
   
    if isstr(dataFile)
        isTif = strsplit(dataFile, '.');
        strIsTif = isTif{2};
        isDcimg =strcmpi(strIsTif, 'dcimg');
    else
        isTif = strsplit(dataFile{1}, '.');
        strIsTif = isTif{2};
        isDcimg =strcmpi(strIsTif, 'dcimg');
    end
            
             avgImg = zeros(imgHeight,imgWidth);
%             avgImgFrames = min(200,length(frames));
            avgImgFrames = min(200,numFrames);
            for a = 1:avgImgFrames
                if isDcimg
%                     dcimgfile = fullfile(dataPath, dataFile_stringname);
%                     dcimgfile = strrep(dcimgfile, '\', '\\');
                    
                    [framedata,totalframes]= dcimgmatlab(a-1, dcimgfile);
                    totalframes = double(totalframes);
                    framedatatrans = transpose (framedata);
                    laserBkgnd =  f_waveletBackground(framedatatrans);
                    %                     avgImg = avgImg + double(framedatatrans) - darkAvg;
%                     avgImg = avgImg + laserBkgnd - darkAvg;
                    avgImg = avgImg + laserBkgnd;

                    
                    
                else
                    dcimgfile = fullfile(dataPath, dataFile);
                    avgImg = avgImg + double(imread([dataFile{stack}],frames(a),'Info', fileInfo)) - darkAvg;
                end
                a
            end
            avgImgLaser = avgImg/avgImgFrames;
            
            threshold = 175;
%             threshold = 160;
            avgImgLaser(avgImgLaser < threshold) = 0;
%             
% %             if ~isequal(size(gain),[length(avgImgLaser(1,:)), length(avgImgLaser(:,1))])
%                 if channel == 'g'
%                     load('gainGreen.mat', 'gain');
%                 elseif channel == 'r'
%                     load('gainRed.mat', 'gain');
%                 else
%                     load('gainRed.mat', 'gain');
%                 end
% %             end
            
            %This is measured in counts
%             bkgndImg_avg = bkgndImgTotal/numbkgndImg;
%               avgImgLaser = avgImgLaser./gain;
           
            %finds laser intensity at each point in the FOV (does not fit
            %gaussian)
            
%             replotLaser = 1;
%             while replotLaser == 1

%             figure;
%             imagesc(avgImgLaser,[0 300]);
%             axis image;colorbar;colormap hot;
%             title('Crop out area of high intensity (Beads)');
%             beadLoc = roipoly;
%             close
%             
%             avgImgLaser(beadLoc == 1) = min(min(avgImgLaser));
            
            %3D laser profile
            laserProfileAvg=avgImgLaser*powerAtObjective/sum(sum(avgImgLaser))...
                /(nmPerPixel^2)*10^14;

%             %2D laser profile                           
%             imagesc(laserProfile);
%             axis square
%             colorbar
%             title('2D Laser Intensity Profile');
            
%             smooth 3D laser profile
            laserFilter = fspecial('disk',21);
            laserProfileSmooth = imfilter(laserProfileAvg,laserFilter);
 
            figure
            surf(laserProfileSmooth)
%             surf(laserProfile)
            shading interp
            title('3D Laser Intensity Profile');
            zlabel('Intensity, W/cm^2');
            
%             prompt = 'Do you want to try to replot laser profile?';
%             def = {'0'};
%             num_lines = 1;
%             dlg_title = 'Replot Laser?';
%             inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
%             replotLaser = str2double(inputdialog(1));     
% %             end
            z = 1;
            
            