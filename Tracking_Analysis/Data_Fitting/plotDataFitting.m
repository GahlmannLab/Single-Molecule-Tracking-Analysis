%plotDataFitting
% This function will plot the PDF of the data as well as the fitted curves

%Load in file saved from Data_Fitting.m
[dataFile1, dataPath1] = uigetfile({'*.mat';'*.*'},'Open file with fitted data','MultiSelect', 'on');
load([dataPath1 dataFile1]);


if isDiff && ~isMinflux
    load('interpCurvesAll4_PDF.mat');
    
    close all;
    xq = 0:0.01:15;
    figure;
    hold on;
    qry = 10;
    xq1 = xq(1:qry:end);
    
    nterm = length(xParamFinal)/2;
    
    a=xParamFinal(1:nterm); % coefficients of basis functions
    u=xParamFinal(nterm+1:end); % parameters for selecting basis functions
    af = x0_Fix(1:ntermf);
    uf = x0_Fix(ntermf+1:end);
    
    sumTerms=zeros(size(xq));
    sumTermsTemp=[];
    
    for k=1:nterm
        temp = [];
        yq = ones(1,length(xq))*u(k);
        
        if nterm == 1
            temp = F(xq,yq)*a(k);
        else
            temp = F(xq,yq)*a(k);
        end
        
        sumTerms = sumTerms + temp;
        temp1 = temp(1:qry:end);
        plot(xq1,temp1);
        
    end
    
    sumTermsFixed = zeros(1,length(sumTerms));
    for j = 1:ntermf
        %add in fixed components
        yq = ones(1,length(xq))*uf(j);
        tempf = F(xq,yq)*af(j);
        
        sumTermsFixed = sumTermsFixed + tempf;
    end
    sumTermsFixed1 = sumTermsFixed(1:qry:end);
    plot(xq1,sumTermsFixed1);
    
    sumTerms = sumTerms + sumTermsFixed;
    

    sumTerms1 = sumTerms(1:qry:end);
    [N,edges] = histcounts(data,xq1);
    plot(xq1,sumTerms1);
    bar((xq1(2:end)-(xq(2)/2)),(N/sum(N))*sum(sumTerms1));
    xlim([0 10]);
    limy=get(gca,'YLim');
    ylim([0 limy(2)]);
    
    annotation('textbox','String',[num2str(xParamFinal) '    ' num2str(x0_Fix_pop)])
    
    annStr = [];
    if x0_Fix_pop ~=0
        annStr{1} = [num2str(x0_Fix_pop) '    <0.5'];
    else
        annStr{1} = [' '];
    end
    for a = 1:nterm
        annStr{a+1} =  [num2str(xParamFinal(a)) '     ' num2str(xParamFinal(nterm+a))];
    end
    t =annotation('textbox','String',annStr,'Position',[0.6, 0.6, 0.1, 0.1]);
    t.FontSize = 14;
    set(gca,'children',flipud(get(gca,'children')))
      
elseif ~isDiff && ~isMinflux
    % Plot Displacements
    load('interpCurvesDisp_PDF.mat');
    
    close all;
    xq = 0:1:1500;
    figure;
    hold on;
    qry = 10;
    xq1 = xq(1:qry:end);
    
    nterm = length(xParamFinal)/2;
    
    a=xParamFinal(1:nterm); % coefficients of basis functions
    u=xParamFinal(nterm+1:end); % parameters for selecting basis functions
    af = x0_Fix(1:ntermf);
    uf = x0_Fix(ntermf+1:end);
    
    sumTerms=zeros(size(xq));
    sumTermsTemp=[];
    
    for k=1:nterm
        % create linear combination of basis functions
        temp = [];
        yq = ones(1,length(xq))*u(k);
        %multiply cdf extracted from scattered interpolant by
        %population fraction parameter. It is possible to add them this
        %way because all cdf's are normalized to 1.
        
        if nterm == 1
            temp = F(xq,yq)*a(k);
        else
            temp = F(xq,yq)*a(k);
        end
        
        sumTerms = sumTerms + temp;
        temp1 = temp(1:qry:end);
        plot(xq1,temp1);
        
    end
    
    sumTermsFixed = zeros(1,length(sumTerms));
    for j = 1:ntermf
        %add in fixed components
        yq = ones(1,length(xq))*uf(j);
        tempf = F(xq,yq)*af(j);
        
        sumTermsFixed = sumTermsFixed + tempf;
    end
    sumTermsFixed1 = sumTermsFixed(1:qry:end);
    plot(xq1,sumTermsFixed1);
    
    sumTerms = sumTerms + sumTermsFixed;
    
    sumTerms1 = sumTerms(1:qry:end);
    [N,edges] = histcounts(data,xq1);
    plot(xq1,sumTerms1);
    bar((xq1(2:end)-(xq(2)/2)),(N/sum(N))*sum(sumTerms1));
    xlim([0 1500]);
    limy=get(gca,'YLim');
    ylim([0 limy(2)]);
        
    annStr = [];
    if x0_Fix_pop ~=0
        annStr{1} = [num2str(x0_Fix_pop) '    <0.5'];
    else
        annStr{1} = [' '];
    end
    for a = 1:nterm
        annStr{a+1} =  [num2str(xParamFinal(a)) '     ' num2str(xParamFinal(nterm+a))];
    end
    t =annotation('textbox','String',annStr,'Position',[0.6, 0.6, 0.1, 0.1]);
    t.FontSize = 14;
    set(gca,'children',flipud(get(gca,'children')))
    
    
elseif ~isDiff && isMinflux
    % Plot Minflux Displacements
    load('interpCurvesMinfluxDisp.mat');
    
    close all;
    xq = 0:1:250;
    figure;
    hold on;
    qry = 5;
    xq1 = xq(1:qry:end);
    
    nterm = length(xParamFinal)/2;
    
    a=xParamFinal(1:nterm); % coefficients of basis functions
    u=xParamFinal(nterm+1:end); % parameters for selecting basis functions
    af = x0_Fix(1:ntermf);
    uf = x0_Fix(ntermf+1:end);
    
    sumTerms=zeros(size(xq));
    sumTermsTemp=[];
    
    for k=1:nterm
        temp = [];
        yq = ones(1,length(xq))*u(k);   
        if nterm == 1
            temp = F(xq,yq)*a(k);
        else
            temp = F(xq,yq)*a(k);
        end
        
        sumTerms = sumTerms + temp;
        temp1 = temp(1:qry:end);
        plot(xq1,temp1);
        
    end
    
    sumTermsFixed = zeros(1,length(sumTerms));
    for j = 1:ntermf
        yq = ones(1,length(xq))*uf(j);
        tempf = F(xq,yq)*af(j);
        
        sumTermsFixed = sumTermsFixed + tempf;
    end
    sumTermsFixed1 = sumTermsFixed(1:qry:end);
    plot(xq1,sumTermsFixed1);
    
    sumTerms = sumTerms + sumTermsFixed;
    
    sumTerms1 = sumTerms(1:qry:end);
    [N,edges] = histcounts(data,xq1);
    plot(xq1,sumTerms1);
    bar((xq1(2:end)-(xq(2)/2)),(N/sum(N))*sum(sumTerms1));
    xlim([0 300]);
    limy=get(gca,'YLim');
    ylim([0 limy(2)]);
        
    annStr = [];
    if x0_Fix_pop ~=0
        annStr{1} = [num2str(x0_Fix_pop) '    <0.5'];
    else
        annStr{1} = [' '];
    end
    for a = 1:nterm
        annStr{a+1} =  [num2str(xParamFinal(a)) '     ' num2str(xParamFinal(nterm+a))];
    end
    t =annotation('textbox','String',annStr,'Position',[0.6, 0.6, 0.1, 0.1]);
    t.FontSize = 14;
    set(gca,'children',flipud(get(gca,'children')))
    
end
