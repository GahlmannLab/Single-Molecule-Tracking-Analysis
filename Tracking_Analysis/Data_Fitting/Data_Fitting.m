% Julian Rocha

%Data Fitting
%This function takes a distribution of diffusion coefficients for molecular
%displacements and fits the distribution with simulated distribution(s).

%Input relevent information.

% Max number of states = highest number of states to attempt to fit. The
% higher this value, the longer the fitting routine will take as it has to
% try more combinations.

% CDF Fitting = 1 if for cumulative distribution function (CDF) or 0 for
% probability density funciton (PDF) fitting. This value should be set to 1
% for a more accurate fitting.

% Diffusion coefficient fitting = 1 for diffusion coefficient distribution
% fitting, or 0 for fitting of the displacement distribution.

% Minflux fitting = 1 for fitting of simulated Minflux displacement
% distribtuions.

input = 'Input';
prompt = {'Enter Max Number of States to Fit','CDF Fitting?','Diffusion Coefficient Fitting?','Minflux fitting?'};
dims = [1 35];
definput = {'3','1','1','0'};
tempAns = inputdlg(prompt,input,dims,definput);

numStatesMax = str2double(tempAns{1});
isCDF = str2double(tempAns{2});
isDiff = str2double(tempAns{3});
isMinflux = str2double(tempAns{4});



% The name of the file must be the same as the variable in the file with the
% diffusion coefficients/displacements
[dataFile1, dataPath1] = uigetfile({'*.mat';'*.*'},'Open file with experimental data','MultiSelect', 'on');

S = load([dataPath1 dataFile1], dataFile1(1:end-4));      % experimental data
data = struct2array(S);

% Remove datapoints > 20 µm^2/s for difussion coefficients, 2500 nm for
% displacement distribtuions, or 250 nm for Minflux data so that curve
% interpolation is monotonically increasing.
if isMinflux
    xMax = 250;
    data = data(data<xMax);
elseif isDiff == 1
    xMax = 20;
    data = data(data<xMax);
elseif ~isDiff
    xMax = 2500;
    data = data(data<xMax);
end

%Load in the appropriate library of simulated distributions based on the
%type of fitting being performed (i.e. diffusion coefficients/displacements
%or CDF/PDF fitting).
if isMinflux
    %this is for CDF fitting of Minflux Displacements
    load('interpCurvesBinlessMinfluxDisp.mat')
    load('designMatrix_MinfluxDisp.mat');
elseif isCDF && isDiff
    load('interpCurvesBinless.mat');
    load('designMatrix_0.05_0_15.mat');
elseif isCDF && ~isDiff
    load('interpCurvesBinless_Disp.mat');
    load('designMatrix_5_5_2000_DISP.mat');
elseif ~isCDF && isDiff
    load('interpCurves_PDF.mat');
    load('designMatrix_0.05_0_15_PDF.mat');
elseif ~isCDF && ~isDiff
    load('interpCurvesDisp_PDF.mat');
    load('designMatrix_5_5_2000_DISP_PDF.mat');
end


% Perform a linear least squares fitting to get an initial guess for the
% number of diffusive states present
[x0_Final x0_Fix numStatesLLS] = LLS_Initial_States(C,dMax,dRange,resolution,xq,data,numStatesMax,isCDF);


%% Peform non-linear least squares fit with k-means cross-validation, using parameter vectors created in LLS_Initial_States

%Calculate the sum of the fixed terms in x0_Fix
sumTermsFix = zeros(size(xq));
ntermf = length(x0_Fix)/2;
af = x0_Fix(1:ntermf);
uf = x0_Fix(ntermf+1:end);
for k = 1:ntermf
    yq = ones(1,length(xq))*uf(k);
    tempf = F(xq,yq)*af(k);
    sumTermsFix = sumTermsFix + tempf;
end
fixPop = sum(x0_Fix(1:ntermf));


if numStatesMax > 1
    %split data into 5 (or 10) random sets for k-means cross-validation
    K = 5;
    N = length(data);
    dataIdx = crossvalind('Kfold', N, K);
    errorCrossVal = cell(1,numStatesMax);
    
    x0_all = x0_Final;
    
    %errorCrossVal contains the error associated with fitting each of the 5
    %(or 10) data sets for each parameter vector in x0_all
    for a = 1:numStatesMax
        errorCrossVal{a} = zeros(length(x0_all{a}),K);
    end
    
    for j = 1:K
        dataTrain = data(dataIdx~=j);
        dataValid = data(dataIdx==j);
        %parfor loop can be used to decrease run-time
        %         parfor a = 1:numStatesMax
        for a = 1:numStatesMax
            for b = 1:length(x0_all{a})
                diffData =  dataTrain;
                x0 = [x0_all{a}{b}(:,2)' x0_all{a}{b}(:,1)'];
                x0(length(x0)/2) = [];
                [xParam errorOut] =  Binless_Fitting_LLS_start(diffData,F,x0,x0_Fix,isCDF,isDiff,isMinflux);
                errorCrossVal{a}(b,j) = errorOut;
            end
        end
    end    
    
    %compute the minumum mean error for each number of states
    errorCrossValMean = cell(numStatesMax,1);
    errorCrossValMin = zeros(1,numStatesMax);
    minIdx = zeros(1,numStatesMax);
    for a = 1:numStatesMax
        for b = 1:length(errorCrossVal{a})
            errorCrossValMean{a} = mean(errorCrossVal{a},2);
            [errorCrossValMin(a) minIdx(a)] = min(errorCrossValMean{a});
        end
    end
    
%     figure;
%     plot(1:numStatesMax,errorCrossValMin);
    
    close all;
    
    %Find point at which minimum does not decrease by at least 5%
    %(errorThresh)
    errorThresh = 0.05;
    for a = 2:numStatesMax
        errorThreshVal = errorCrossValMin(a-1) - errorThresh*errorCrossValMin(a-1);
        if errorCrossValMin(a) > errorThreshVal
            minNumStates = a-1;
            break;
        end
    end
    
    if ~exist('minNumStates')
        minNumStates = numStatesMax;
    end
    %Create the final fitting parameter vector x0_Final by finding the
    %parameter set with the lowest error in the errorCrossValMin variable
    x0_Final = [x0_all{minNumStates}{minIdx(minNumStates)}(:,2)' x0_all{minNumStates}{minIdx(minNumStates)}(:,1)'];
    x0_Final(length(x0_Final)/2) = [];
else
    x0_Final = x0_Final{1}{1};
    x0_Final(length(x0_Final)) = [];
end

%Use the final fitting parameter vector x0_Final to do a non-linear fitting
%on the full data set.
[xParamFinal errorOutFinal] =  Binless_Fitting_LLS_start(data,F,x0_Final,x0_Fix,isCDF,isDiff,isMinflux)



x0_Fix_pop = sum(x0_Fix(1:length(x0_Fix)/2));
if isDiff == 1
    save([dataPath1 date '_DIFF_FinalFit_' dataFile1]);
else
    save([dataPath1 date '_DISP_FinalFit_' dataFile1]);
end



%% Bootstrap Results

msRun = 1;
nboot = 2;
bootstatTot = [];
for a = 1:50
%Give a slight offset to the starting parameters to avoid getting stuck in
%a local minima when fitting with the same starting parameters
x0_Final2 = (x0_Final+x0_Final*0.2.*(rand(1,3)-0.5)*2);
[bootstat, bootsam] = bootstrp(nboot,@(x)Binless_Fitting_LLS_start(x,F,x0_Final2,x0_Fix,isCDF,isDiff,isMinflux),data);
end
bootstatTot = vertcat(bootstatTot,bootstat);
meanParam = mean(bootstatTot);
stdParam = std(bootstatTot);

if isDiff == 1
    save([dataPath1 date '_DIFF_bootstrap_' dataFile1]);
else
    save([dataPath1 date '_DISP_bootstrap_' dataFile1]);
end


