
function [x0_Final x0_Fix numStatesLLS] = LLS_NumStates_Initial(C,dMax,dRange,resolution,xq,data,numStatesMax,isCDF)
% This fucntions perform a linear least squares fitting to get an initial guess for the
% number of diffusive states present

data = data;
%% Perform an interpolation of the data
if isCDF == 1
[curveInterp xData data_interp] = interpCurve(data,xq);
else
[curveInterp xData data_interp] = interpCurvePDF(data,xq);
data_interp = data_interp/sum(data_interp);
end

%% Solve the constrained linear least square problem with the full dataset
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter-detailed', ...
    'MaxIter',1000, 'TolFun', 1e-15);

Aeq = ones(1,length(dRange));
beq = 1;

lb = zeros(length(dRange),1);
ub = ones(length(dRange),1);

% tic
linearCoeff_full = lsqlin(C,data_interp',[],[],Aeq,beq,lb,ub,[],options);
% toc


out = dRange(find(linearCoeff_full>0.005))';
out = [out , linearCoeff_full(linearCoeff_full>0.005)];


%% Combine and separate for different number of state fits

%first step is removing any peak that is less than 5% of the total
%population and adding it to the closest peak. This is done to simplify the
%fitting.
outFix = out;
minThresh = 0.05;
if sum(outFix(:,2) < minThresh)
    fin = 0;
    loopNum = 0;
    while fin == 0
        for a = 1:length(outFix(:,2))
            finTemp = 0;
            if outFix(a,2) < minThresh
                if outFix(a,1) <= 0.5
                    b = 1;
                    if a == 1
                        if outFix(a+b,1) < 0.5
                            outFix(a+b,2) = outFix(a+b,2) + outFix(a,2);
                            outFix(a+b,1) = outFix(a,1)*(outFix(a,2)/(outFix(a+b,2)+outFix(a,2)))...
                                + outFix(a+b,1)*(outFix(a+b,2)/(outFix(a+b,2)+outFix(a,2)));
                            outFix(a,:) = [];
                            finTemp = 1;      
                        else
                            finTemp = 1;
                        end
                    else
                        if outFix(a+b,1) < 0.5
                            outFix(a+b,2) = outFix(a+b,2) + outFix(a,2);
                            outFix(a+b,1) = outFix(a,1)*(outFix(a,2)/(outFix(a+b,2)+outFix(a,2)))...
                                + outFix(a+b,1)*(outFix(a+b,2)/(outFix(a+b,2)+outFix(a,2)));
                            outFix(a,:) = [];
                            finTemp = 1; 
                        elseif outFix(a-b,1) < 0.5
                            outFix(a-b,2) = outFix(a-b,2) + outFix(a,2);
                            outFix(a-b,1) = outFix(a,1)*(outFix(a,2)/(outFix(a-b,2)+outFix(a,2)))...
                                + outFix(a-b,1)*(outFix(a-b,2)/(outFix(a-b,2)+outFix(a,2)));
                            outFix(a,:) = [];
                            finTemp = 1; 
                        else
                            finTemp = 1;
                        end
                    end
                else
                    if a == 1
                        b = 1;
                        while finTemp == 0
                            if outFix(a+b,2) >= minThresh
                                outFix(a+b,2) = outFix(a+b,2) + outFix(a,2);
                                outFix(a,:) = [];
                                finTemp = 1;
                            else
                                b = b+1;
                            end
                        end
                    elseif a == length(outFix(:,2))
                        b = 1;
                        while finTemp == 0
                            if outFix(a-b,2) >= minThresh
                                outFix(a-b,2) = outFix(a-b,2) + outFix(a,2);
                                outFix(a,:) = [];
                                finTemp = 1;
                                fin = 1;
                            else
                                b = b+1;
                            end
                        end
                    else
                        b = 1;
                        while finTemp == 0
                            if outFix(a-b,2) >= minThresh && outFix(a+b,2) >= minThresh
                                diff1 = outFix(a,1)-outFix(a-b,1);
                                diff2 = outFix(a+b,1)-outFix(a,1);
                                if diff1 < diff2
                                    outFix(a-b,2) = outFix(a-b,2) + outFix(a,2);
                                    outFix(a,:) = [];
                                    finTemp = 1;
                                else
                                    outFix(a+b,2) = outFix(a+b,2) + outFix(a,2);
                                    outFix(a,:) = [];
                                    finTemp = 1;
                                end
                            elseif outFix(a-b,2) >= minThresh
                                outFix(a-b,2) = outFix(a-b,2) + outFix(a,2);
                                outFix(a,:) = [];
                                finTemp = 1;
                            elseif outFix(a+b,2) >= minThresh
                                outFix(a+b,2) = outFix(a+b,2) + outFix(a,2);
                                outFix(a,:) = [];
                                finTemp = 1;
                            else
                                b = b+1;
                            end
                        end
                    end
                end
            elseif a == length(outFix(:,2)) && outFix(a,2) >= minThresh
                fin = 1;
            end
            loopNum = loopNum+1;
            if loopNum >= 100
                fin = 1;
            end
            if finTemp == 1
                break;
            end
        end
    end
end

%Next, if there are components with d* < 0.50 µm^2/s, hold them constant
outFix2 = outFix(outFix(:,1)>0.5,:);
x0_Fix_temp = outFix(outFix(:,1)<=0.5,:);
x0_Fix = [x0_Fix_temp(:,2)' x0_Fix_temp(:,1)'];


%Next, combine diffusive states or separate them to get parameter vectors
%with different lengths
xParamFull = cell(1,numStatesMax);
errorOutFull = zeros(1,numStatesMax);
x0_mod_tot = cell(1,numStatesMax);
x0 = outFix2;
numStatesOrig = (length(x0))/2;
percThresh = 0.2;


%create initial parameter vector by combining LLS states within 20%
%(percThresh variable). Again this is done to simplify the initial fitting.
x0_temp = x0;
tempIdx = zeros(length(x0_temp(:,1)),1);
%tempIdx, 1=state above is within percThresh, 2=state below is within
%percThresh
for c = 1:length(x0_temp(:,1))
    if c ==1
        if length(x0_temp(:,1)) == 1
            break;
        elseif x0_temp(c+1,1) - x0_temp(c,1) <= (x0_temp(c,1)*percThresh)
            tempIdx(c) = 1;
            tempIdx(c+1) = 2;
        end
    elseif c == length(x0_temp)
        if x0_temp(c,1) - x0_temp(c-1,1) <= (x0_temp(c,1)*percThresh)
            tempIdx(c-1) = 1;
            tempIdx(c) = 2;
        end
    else
        if x0_temp(c,1) - x0_temp(c-1,1) <= (x0_temp(c,1)*percThresh)
            tempIdx(c) = 2;
            tempIdx(c-1) = 1;
        end
        if x0_temp(c+1,1) - x0_temp(c,1) <= (x0_temp(c,1)*percThresh)
            tempIdx(c) = 1;
            tempIdx(c+1) = 2;
        end
    end
end

if sum(tempIdx) > 0
    fin = 0;
    while fin == 0
        for c = 1:length(x0_temp)
            if tempIdx(c) == 1
                x0_temp(c+1,2) = x0_temp(c,2) + x0_temp(c+1,2);
                x0_temp(c+1,1) = sum(x0_temp(c:c+1,1).*(x0_temp(c:c+1,2)/sum(x0_temp(c:c+1,2))));
                x0_temp(c,:) = -1*ones(1,length(x0_temp(c,:)));
                tempIdx(c) = -1;
                if c == length(x0_temp)
                    fin = 1;
                end
                break;
            elseif tempIdx(c) ~= 1 && c == length(x0_temp)
                fin = 1;
            end
            
        end
    end
end
%x0_temp is the weighted average for components with d within 20% of each
%other
tempIdx = tempIdx(tempIdx >= 0);
x0_temp = x0_temp(x0_temp(:,1) >= 0,:);

%Combine or separate states in the inital parameter vector to get different
%combinations of states for each numbers of states
numStatesLLS = length(x0_temp(:,1));
numLower = numStatesLLS-1;
for a = numLower:-1:1
    if a == 1
        temp = x0_temp;
        temp(1,1) = sum(x0_temp(:,1).*(x0_temp(:,2)/sum(x0_temp(:,2))));
        temp(1,2) = sum(x0_temp(:,2));
        x0_mod_tot{1}{1} = temp(1,:);
    elseif a == numLower
        for b = 1:a
            temp = [];
            temp = x0_temp;
            temp(b+1,1) = sum(x0_temp(b:b+1,1).*(x0_temp(b:b+1,2)/sum(x0_temp(b:b+1,2))));
            temp(b+1,2) = sum(x0_temp(b:b+1,2));
            temp(b,:) = [];
            x0_mod_tot{a}{b} = temp;
        end
        
    else
        d = 1;
        for c = 1:length(x0_mod_tot{a+1})
            for b = 1:a
                temp = [];
                temp = x0_mod_tot{a+1}{c};
                temp2 = x0_mod_tot{a+1}{c};
                temp(b+1,1) = sum(temp2(b:b+1,1).*(temp2(b:b+1,2)/sum(temp2(b:b+1,2))));
                temp(b+1,2) = sum(temp2(b:b+1,2));
                temp(b,:) = [];
                x0_mod_tot{a}{d} = temp;
                d = d+1;
            end
        end
    end
    
end

%remove duplicates in x0_mod_tot
for a = 2:numLower-1
    for b = 1:length(x0_mod_tot{a})-1
        for c = b+1:length(x0_mod_tot{a})
            if ~isempty(x0_mod_tot{a}{b}) && ~isempty(x0_mod_tot{a}{c})
                if sum(sum(round(x0_mod_tot{a}{b},3) == round(x0_mod_tot{a}{c},3))) == length(x0_mod_tot{a}{b})*2
                    x0_mod_tot{a}{c} = [];
                end
            end
        end
    end
    x0_mod_tot{a} = x0_mod_tot{a}(~cellfun('isempty',x0_mod_tot{a})) ;
end
x0_mod_tot{numStatesLLS}{1} = x0_temp;


for a = numStatesLLS+1:numStatesMax
    if a == numStatesLLS+1
        x0_mod_tot{a} = cell(a-1,1);
        for b = 1:a-1
            temp = zeros(length(x0_temp(:,1))+1,length(x0_temp(1,:)));
            temp(1:b,:) = x0_temp(1:b,:);
            temp(b+2:end,:) = x0_temp(b+1:end,:);
            temp(b,1) = x0_temp(b,1) - x0_temp(b,1)*0.2;
            temp(b+1,1) = x0_temp(b,1) + x0_temp(b,1)*0.2;
            temp(b:b+1,2) = x0_temp(b,2)/2;
            x0_mod_tot{a}{b} = sortrows(temp,1);
        end
    else
        d = 1;
        for b = 1:a-2
            x0_temp2 = x0_mod_tot{a-1}{b};
            for c = 1:a-1
                temp = zeros(length(x0_temp2(:,1))+1,length(x0_temp2(1,:)));
                temp(1:c,:) = x0_temp2(1:c,:);
                temp(c+2:end,:) = x0_temp2(c+1:end,:);
                temp(c,1) = x0_temp2(c,1) - x0_temp2(c,1)*0.2;
                temp(c+1,1) = x0_temp2(c,1) + x0_temp2(c,1)*0.2;
                temp(c:c+1,2) = x0_temp2(c,2)/2;
                x0_mod_tot{a}{d} = sortrows(temp,1);
                d = d + 1;
            end
        end
    end
    
    
end

for a = numStatesLLS+1:numStatesMax
    %Remove parameter vectors that have split a state into states with
    %population fractions less than the threshold in minThresh
    for e = 1:length(x0_mod_tot{a})
        if sum(x0_mod_tot{a}{e}(:,2) < minThresh) ~= 0
            x0_mod_tot{a}{e} = [];
        end
    end
    x0_mod_tot{a} = x0_mod_tot{a}(~cellfun('isempty',x0_mod_tot{a}));
end

%x0_Final is the final array of fitting parameter vectors. It is a cell
%array with a length equal to numStatesMax. Each cell contains a cell array
%of different parameter vectors with a given number of states.
x0_Final = x0_mod_tot;

if isempty(x0_Fix)
    x0_Fix = [0 0];
end
end

