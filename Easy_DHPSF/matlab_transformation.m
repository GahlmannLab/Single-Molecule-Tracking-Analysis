function [tform, FRE, TRE, FRE_full, TRE_full] = matlab_transformation(...
    cp_channel1, cp_channel2 , tform_mode)
% This function computes a 2D transformation using the given set of control points
% and also computes the associated FRE and TRE values.
% The code is adapted from "Single-Molecule High Resolution Colocalization of Single
% Probes" by L. Stirling Churchman and James A. Spudich in Cold Spring
% Harbor Protocols (2012), doi: 10.1101/pdb.prot067926
% Modified Definitions of FRE and TRE - Andreas Gahlmann, 20110511

% Transform the data using the cp2tform command
tform = cp2tform(cp_channel1,cp_channel2,tform_mode);

% Calculate the metrices to estimate the error associated with this
% Calculate the fiducial registration error (FRE)
trans_cp_channel2 = tforminv(cp_channel2,tform);

% FRE = sqrt(sum(sum((cp_channel1-trans_cp_channel2).^2))/(length(cp_channel1)));
% FRE_full = ((cp_channel1-trans_cp_channel2).^2)/(length(cp_channel1));

FRE_full = sqrt(sum((cp_channel1-trans_cp_channel2).^2,2));
FRE = mean(FRE_full)

% Calculate the target registration error (TRE)
number_cp = length(cp_channel1); % find the number of control points
% Loop through the control points
for i=1:number_cp
    i
    remain_cp = [1:i-1 i+1:number_cp]; % take out that control point
    
    % Calculate the transformation without the ith control point
    tform = cp2tform(cp_channel1(remain_cp,:),cp_channel2(remain_cp,:),tform_mode);
    
    % Transform left out control point with above found transformation
    trans_cp_channel2(i,:) = tforminv(cp_channel2(i,:),tform);
    
end

% TRE = sqrt(sum(sum((cp_channel1 - trans_cp_channel2).^2))/(length(cp_channel1)));
% TRE_full = ((cp_channel1 - trans_cp_channel2).^2)/(length(cp_channel1));
TRE_full = sqrt(sum((cp_channel1-trans_cp_channel2).^2,2));
TRE = mean(TRE_full(TRE_full<100000));

% Restore the full transform function again
tform = cp2tform(cp_channel1,cp_channel2,tform_mode);
channel2_trans = tforminv(cp_channel2,tform);

% show the results
figure
distlimit = 5;

subplot(2,2,1)
scatter(cp_channel1(:,1), cp_channel1(:,2))
title({'Reflected Channel';'Channel 1'})
hold on
scatter(channel2_trans(:,1), channel2_trans(:,2), 10, 'filled')
hold off
axis square
subplot(2,2,2)
scatter(cp_channel2(:,1), cp_channel2(:,2))
title({'Transmitted Channel';'Channel 2'})
axis square

subplot(2,2,3)
hist(FRE_full(FRE_full<=distlimit), 20)
title({['Target Registration Error']; [tform_mode ' Transformation']});
xlabel('Distance (nm)');
ylabel('Frequency');
xlim([0 distlimit]);
legend(['Mean = ' num2str(FRE, 3) ' pixel']);

subplot(2,2,4)
hist(TRE_full(TRE_full<distlimit),20)
title({['Fiducial Registration Error']; [tform_mode ' Transformation']});
xlabel('Distance (nm)');
ylabel('Frequency');
xlim([0 distlimit]);
legend(['Mean = ' num2str(TRE, 3) ' pixel']);

% saveas(gcf,['Tranformation_' tform_mode '.fig']);
% saveas(gcf,['z0nm_stats_' tform_mode '.png']);

end