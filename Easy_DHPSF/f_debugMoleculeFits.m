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

function f_debugMoleculeFits(totalPSFfits,numFrames,dataFile,dataPath,...
            darkFile,darkPath,ROI,templateFrames,peakThreshold,threshFile,...
            gaussianFilterSigma,minDistBetweenSMs)
    %% initialize variables, output path, etc.

    scrsz = get(0,'ScreenSize');
    
    if mod(ROI(3),2)==1
        ROI(3) = ROI(3)-1;
    end
    if mod(ROI(4),2)==1
        ROI(4) = ROI(4)-1;
    end
    cropWidth = ROI(3);
    cropHeight = ROI(4);
    
    [debugFile, debugPath] = uiputfile({'*.*'},'Specify a path and folder name for output');
    if isequal(debugFile,0)
        return;
    end
    outputPrefix = [debugPath filesep debugFile filesep];
    outputPrefixImages = [outputPrefix filesep 'images' filesep];
    mkdir(outputPrefix);
    mkdir(outputPrefixImages);
    
    %% deal with templates, to get template matches in the output
    
    if exist('threshFile')
    numPSFfits = 0;
    load(threshFile,'template');
    templateSize = size(template,2);
    numTemplates = length(templateFrames);
    templateColors = jet(numTemplates);
    
    templatePad = zeros(numTemplates,cropHeight,cropWidth);
    templateFT = zeros(numTemplates,cropHeight,cropWidth);
        for a=1:numTemplates
%             
%             if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
%                 templatePad(a,:,:) = padarray(squeeze(template(a,:,:)),...
%                     [(cropHeight-size(template,2))/2 ...
%                     (cropWidth-size(template,3))/2],min(min(template(a,:,:))));
%             else
                templatePad(a,:,:) = padarray(squeeze(template(templateFrames(a),:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(templateFrames(a),:,:))));
%             end
            
            % multiplying by conjugate of template in FT domain is equivalent
            % to flipping the template in the real domain
            templateFT(a,:,:) = conj(fft2(squeeze(templatePad(a,:,:))));
        end
        clear templatePad;
        
        % apply Gaussian filter to phase correlation data to weight low frequencies
        % more heavily since SNR is higher there
        gaussianFilter = abs(fft2(fspecial('gaussian', [cropHeight cropWidth], gaussianFilterSigma)));
    end    
    
    %% spits out selected images to debug whether fitting appropriate

    % create dark offset array ('darkAvg')
    if ~isequal(darkFile,0)
        % Computes average of darkAvg frames for background subtraction
        darkFileInfo = imfinfo(darkFile);
        numDarkFrames = length(darkFileInfo);
        darkAvg = zeros(darkFileInfo(1).Height,darkFileInfo(1).Width);
        for dframe = 1:numDarkFrames
            darkAvg = darkAvg + double(imread(darkFile,dframe,'Info',darkFileInfo));
        end
        darkAvg = darkAvg/numDarkFrames;
%         if ~isequal(size(darkAvg),[imgHeight imgWidth])
%             warning('Dark count image and data image stack are not the same size. Resizing dark count image...');
%             darkAvg = imresize(darkAvg,[imgHeight imgWidth]);
%         end
    else
        darkAvg = 0;
    end
    clear darkFileInfo;
    
    % choose frames to output
    startFrame = min(totalPSFfits(:,1)); % this includes fits from all files
    endFrame = max(totalPSFfits(:,1));

    dlg_title = 'Please choose frames to visualize';
    prompt = {  'Frames:'...
                'Print output into current directory?'...
    };
    def = { ...
    num2str(startFrame:100:endFrame), ...
    'Yes'
    };
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    % fitParam = totalPSFfits(:,1:8);
    frames=str2num(inputdialog{1});
    printOutputFrames = strcmp(inputdialog{2},'Yes');
    
    % split chosen frames up into appropriate files and output
    fileFrames = cumsum(numFrames); % numFrames lists starting frame of each file
    
    for i = 1:length(fileFrames)
    fileInfoAll{i} = imfinfo([dataPath dataFile{i}]); % for all files
    end
    
    % actually go through and output the frames with reconstructions
    
    hSMFits=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    
    for i = frames % frames selected by user
        for j = 1:length(fileFrames)
            if i < fileFrames(j)
                fileNum = j; % now we know the right file
                break % out of the 'j' for loop (already found right file)
            end
        end
        fileInfo = fileInfoAll{fileNum};
        if fileNum > 1
            frame = i - fileFrames(fileNum-1); % frame # within the file
        else
            frame = i; % first file does not need to be corrected
        end
        
        %% generate data (assume wavelet corrected)
        
        waveletCorr = true;
        
        data = double(imread([dataPath dataFile{fileNum}],frame,'Info',fileInfo))-darkAvg;
        data = data(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
        
        if waveletCorr
            bkgndImg_curr = f_waveletBackground(data); % assume you are using the wavelet subtraction
            dataWCorr = data - bkgndImg_curr;
        else
            bkgndImg = zeros(length(ROI(2):ROI(2)+ROI(4)-1),...
                length(ROI(1):ROI(1)+ROI(3)-1));
        end
        %% do template matching
        
        if exist('threshFile')

            if waveletCorr
                dataFT = fft2(dataWCorr,cropHeight,cropWidth);
            else
                dataFT = fft2(data,cropHeight,cropWidth);
            end
            maxPeakImg = zeros(cropHeight,cropWidth);
            % matrix PSFLocs stores information about double helices that were
            % found via template matching
            % rows are different matches
            % [xLocation yLocation matchingTemplateNumber matchConfidence];
            PSFLocs = zeros(100,4);
            numPSFLocs = 0;
            for b=1:numTemplates

                H = gaussianFilter./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));

                peakImg = ifftshift(ifft2(dataFT.*squeeze(templateFT(b,:,:)).*H));

                % normalize response of peakImg by dividing by number of pixels in
                % data
                %peakImg = peakImg / (cropHeight*cropWidth);
                maxPeakImg = max(maxPeakImg, peakImg);

                %threshold = mean(peakImg(:))+peakThreshold*std(peakImg(:));
                peakImg(peakImg < peakThreshold(fileNum,b)) = peakThreshold(fileNum,b);

                if isreal(peakImg) && sum(sum(isnan(peakImg)))==0
                    temp = find(imregionalmax(peakImg));
                else
                    warning('got nans in template match!');
                    peakImg(isnan(peakImg)) = 0; %inserted to deal with NaNs -AC 6/22
                    temp = find(imregionalmax(real(peakImg)));
                end


                % make sure threshold didn't eliminate all peaks and create
                % lots of matches
                if length(temp) < cropHeight*cropWidth/2;
                    [tempY, tempX] = ind2sub([cropHeight cropWidth],temp);
                    PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                        [tempX tempY b*ones(length(temp),1) peakImg(temp)];
                    numPSFLocs = numPSFLocs+length(temp);
                end
            end
            clear H dataFT peakImg

            %% filter out extraneous matches due to very strong signals

            if numPSFLocs > 0
                % sort location matrix in decending order of confidence
                temp = sortrows(PSFLocs(1:numPSFLocs,:),-4);
                % copy most confident match to list of locations
                PSFLocs(1,:) = temp(1,:);
                numPSFLocs = 1;
                for b=2:size(temp,1)
                    % make sure that this candidate location is a minimum distance away
                    % from all other candidate locations
                    if sum((temp(b,1)-PSFLocs(1:numPSFLocs,1)).^2 + (temp(b,2)-PSFLocs(1:numPSFLocs,2)).^2 >= minDistBetweenSMs^2) == numPSFLocs
                        % add it to list of locations
                        numPSFLocs = numPSFLocs + 1;
                        PSFLocs(numPSFLocs,:) = temp(b,:);
                    end
                end
            end
        end
        % generate reconstruction from totalPSFfits

        reconstructImg = bkgndImg_curr;%zeros(cropHeight, cropWidth);
        
        [xIdx, yIdx] = meshgrid(1:cropWidth,1:cropHeight);
        molsInFrame = (totalPSFfits(:,1) == i);
        for j = find(molsInFrame)'
            fitParam = totalPSFfits(j,7:14); % see f_fitSMs
            % shift coordinates to be relative to ROI, not entire dataset
            fitParam(3) = fitParam(3) - (ROI(1)-1);
            fitParam(4) = fitParam(4) - (ROI(2)-1);
            fitParam(5) = fitParam(5) - (ROI(1)-1);
            fitParam(6) = fitParam(6) - (ROI(2)-1);
            reconstructImg = reconstructImg + ...
                fitParam(1).*exp( -((xIdx-fitParam(3)).^2+(yIdx-fitParam(4)).^2.) / (2.*fitParam(7).^2)) ...
                +fitParam(2).*exp( -((xIdx-fitParam(5)).^2+(yIdx-fitParam(6)).^2.) / (2.*fitParam(8).^2));
        end
        
        %% display the data and reconstruction
        set(0,'CurrentFigure',hSMFits);
        
        if exist('threshFile')
            subplot('Position',[0.025 0.025 .85/3 .95]);
            imagesc(maxPeakImg,[0 3*min(peakThreshold(fileNum,:))]);axis image;
            title({'Peaks correspond to likely template matches' ...
                [num2str(numPSFLocs) ' matches found']});
        end
        
        if exist('threshFile')
            subplot('Position',[0.075+.85/3 0.025 .85/3 .95]);
        else
            subplot(1,2,1)
        end
        imagesc(data);axis image;colormap hot;
        if exist('threshFile')
            hold on;
            for b=1:numPSFLocs
                %plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
                %    'MarkerSize', 15*PSFLocs(b,4)/peakThreshold(b), ...
                %    'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
                plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
                    'MarkerSize', 15*PSFLocs(b,4)/peakThreshold(fileNum,PSFLocs(b,3)), ...
                    'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
            end
            hold off;
        end
        title({['Frame ' num2str(i) ': raw data - darkAvg counts'] ...
            ['ROI [xmin ymin width height] = ' mat2str(ROI)]});
        
        if exist('threshFile');
            subplot('Position',[0.125+2*.85/3 0.025 .85/3 .95]);
        else
            subplot(1,2,2)
        end
        %         imagesc(reconstructImg+bkgndMean,[min(data(:)) max(data(:))]);axis image;
        imagesc(reconstructImg,[min(data(:)) max(data(:))]);axis image;
        title({'Image reconstructed from fitted matches' ...
            [num2str(sum(molsInFrame)) ' successful fits']});
        
        drawnow;
        
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
            'Callback','uiresume(gcbf)');
        uiwait(gcf);
        
        if printOutputFrames == 1
            set(gcf,'PaperPositionMode','auto');
            saveas(hSMFits, [outputPrefixImages 'frame ' num2str(frame) '.tif']);
        end
        
        clear fileNum
    end
    close(hSMfits)
    %% generates a .csv and/or .mat containing fits, + a histogram of errors
    
    
    
    edges = [-1007:-1000,-3,1,inf];
    vector = histc(totalPSFfits(:,17),edges);

    hErrors=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');

    bar(vector(1:end-1))
    set(gca,'XTickLabel',{'LS error', 'amp ratio', 'lobe dist',...
                          'sig ratio', 'sig size', 'pks out',...
                          'amp< 0', 'guess out','fit err',...
                          'good fit'})
    set(gca,'FontSize',8);

%     print(hErrors,'-dpng',[outputPrefixImages 'outcomes 1.png']);
    saveas(hErrors,[outputPrefixImages 'outcomes.png']);
    
    %% print .csv
    % open a file for writing
%     [fid,message] = fopen([debugPath debugFile], 'w');
%     if ~isempty(message)
%         error([debugPath debugFile ': ' message]);
%         %return;
%     end
%     % print a title, followed by a blank line
%     fprintf(fid, ['frame num,fit flag,SM idx in frame,template x (pix),template y (pix), template idx,' ...
%         'match strength,amp1,amp2,peak1 x (pix),peak1 y (pix),' ...
%         'peak2 x (pix),peak2 y (pix), sigma1 (pix),sigma2 (pix),mean bkgnd photons,'...
%         'fit error,molecule x (nm),molecule y (nm),DHPSF angle,' ...
%         'num photons,interlobe distance,amplitude ratio,sigma ratio,x fid-corrected (nm),y fid-corrected (nm), z fid-corrected (nm),'...
%         'photons detected,mean background photons\n']);
%     fclose(fid);
% 
%     dlmwrite([debugPath debugFile],totalPSFfits(:,[1 17 2:16 18:end]),'-append');
    save([outputPrefix  'debug output.mat']);
    disp('Debug files written successfully.');
end