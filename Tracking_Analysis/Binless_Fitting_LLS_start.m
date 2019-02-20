function [xParam errorOut] =  Binless_Fitting_LLS_start(diffData,F,x0,x0_Fix,isCDF,isDiff,isMinflux)

% Binless_Fitting uses a linear combination of basis functions from a set of
% basis functions to to curve fit diffData.
% Inputs
%    diffData: diffusion coefficient values
%    F: 2D scattered interpolant of simulated diffusion coefficient curves used for fitting
%    x0: parameter vector for linear combinations of basis functions.
%        Each term in the linear combination has two parameters: a
%        multiplicative coefficient and an index of the basis function used
%        for that term.(ex 1 term: x0=[1 1.2] 2 term: x0=[0.5 0.5 1.2
%        2.5]
%
% Outputs
%    xParam: final output parameter vector
%    errorOut: parameter minimized during fitting optimization

%% Set up initial parameters and interpolate original data with B-spline

%xq is the query points vector in x
if isMinflux
    xq = 0:1:250;
    maxX = 250;
elseif isDiff
    xq = 0:0.05:20;
    maxX = 20;
else
    xq = 0:5:2500;
    maxX = 2500;
end

%interpolate data
if isCDF == 1
    [curveInterp xData data_interp] = interpCurve(diffData,xq);
else
    [curveInterp xData data_interp] = interpCurvePDF(diffData,xq);
    data_interp = data_interp/sum(data_interp);
end

% %query the gridded interpolant curveInterp at xq
xData = xq;
if isCDF
    yData = curveInterp(xq);
else
    yData = curveInterp(xq)/sum(curveInterp(xq));
end

%% Set up fitting optimization

%nterm is the number of species being fit to the curve
nterm = (length(x0)+1)/2;
%ntermf is the number of fixed species
ntermf = length(x0_Fix)/2;

%set lower bounds of diffusion coefficients to 0.05 µm^2/s because we can
%not distinguish any value lower than that

%set lower bounds of population fractions to 0.05 (5%) again because we can
%not reliably distinguish populations low in abundance

%setting upper bounds helps to speed up optimization
lb = [zeros(1,nterm-1) (x0(nterm:end) - x0(nterm:end)*0.3)];
ub = [ones(1,nterm-1)*(1-sum(x0_Fix(1:ntermf))) min((x0(nterm:end)+x0(nterm:end)*0.3),15)];

% Set optimization parameters
options = optimoptions('fmincon', ...
    'Display','off','TolX',1e-18, ...
    'TolFun',1e-18,'MaxFunEvals',10000,'MaxIter',1000);
problem = createOptimProblem('fmincon','objective',@(x)linearComb(x),'x0',x0,...
    'lb',lb,'ub',ub,'options',options);

[xParam fval] = fmincon(problem);

    function errorOut=linearComb(x)
        % Objective function linearComb takes n coefficient parameters and n peak location
        % parameters for n terms in the linear combination
        if nterm == 1
            a = 1;
            u = x(1);
            if exist('x0_Fix')
                af = x0_Fix(1:ntermf);
                uf = x0_Fix(ntermf+1:end);
            end
        else
            a=x(1:nterm-1); % coefficients of basis functions
            u=x(nterm:end); % parameters for selecting basis functions
            af = x0_Fix(1:ntermf);
            uf = x0_Fix(ntermf+1:end);
        end
        sumTerms=zeros(size(xq));
        sumTermsTemp=[];
        
        if isCDF == 1
            for k=1:nterm
                % create linear combination of basis functions
                temp = [];
                yq = ones(1,length(xq))*u(k);
                %multiply cdf extracted from scattered interpolant by
                %population fraction parameter. It is possible to add them this
                %way because all cdf's are normalized to 1.
                
                if nterm == 1
                    temp = F(xq,yq)*a(k);
                elseif k == nterm
                    temp = F(xq,yq)*(1-sum(a(1:nterm-1)) - sum(af(1:ntermf)));
                else
                    temp = F(xq,yq)*a(k);
                end
                
                sumTerms = sumTerms + temp;
                
            end
            
            for j = 1:ntermf
                %add in fixed components
                yq = ones(1,length(xq))*uf(j);
                tempf = F(xq,yq)*af(j);
                
                sumTerms = sumTerms + tempf;
            end
        else
            for k=1:nterm
                % create linear combination of basis functions
                temp = [];
                yq = ones(1,length(xq))*u(k);
                %multiply cdf extracted from scattered interpolant by
                %population fraction parameter. It is possible to add them this
                %way because all cdf's are normalized to 1.
                
                if nterm == 1
                    temp = (F(xq,yq)*a(k))/F(xq,yq);
                elseif k == nterm
                    temp = (F(xq,yq)*(1-sum(a(1:nterm-1)) - sum(af(1:ntermf))))/F(xq,yq);
                else
                    temp = (F(xq,yq)*a(k))/F(xq,yq);
                end
                
                sumTerms = sumTerms + temp;
                
            end
            
            for j = 1:ntermf
                %add in fixed components
                yq = ones(1,length(xq))*uf(j);
                tempf = (F(xq,yq)*af(j))*F(xq,yq);
                
                sumTerms = sumTerms + tempf;
            end
        end
        
        %errorOut is residual sum of squares
        errorOut=sum((sumTerms-yData).^2); 
    end

xParam_temp = zeros(1,length(xParam)+1);
xParam_temp(1,nterm) = 1 - sum(xParam(1,1:nterm-1)) - sum(x0_Fix(1:ntermf));
xParam_temp(1,1:nterm-1) = xParam(1,1:nterm-1);
xParam_temp(1,nterm+1:end) = xParam(1,nterm:end);

xParam = xParam_temp;
if isempty(fval)
    fval = 0;
end
errorOut = fval;

end


