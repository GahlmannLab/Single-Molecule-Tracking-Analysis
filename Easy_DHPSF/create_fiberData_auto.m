function dataFileFiber = create_fiberData_auto(savePath)
%     [dataFile, dataPath] = uigetfile({'*.mat';'*.*'},'Open file with localizations');
    
    dataFile = 'Output.mat';
    dataPath = savePath;
    load ([dataPath,dataFile]);

    

    a = 1;
    fiberData(a).xLoc = xLoc;
    fiberData(a).yLoc = yLoc;
    fiberData(a).zLoc = zLoc_IndexCorrected;
    fiberData(a).sigmaX = sigmaX;
    fiberData(a).sigmaY = sigmaY;
    fiberData(a).sigmaZ = sigmaZ;
    fiberData(a).numPhotons = numPhotons;
    fiberData(a).meanBkgnd = meanBkgnd;
    fiberData(a).Frame = frameNum;
%     fiberData1(a).meshPtsX = meshPts_trans2{a,1}(:,1);
%     fiberData1(a).meshPtsY = meshPts_trans2{a,1}(:,2);
    fiberData(1).FrameRange = frameRange;
    
    dVal = strsplit(dataPath,'d');
    dVal = dVal{2};
    dVal = strsplit(dVal,'_');
    dVal = dVal{1};
    fiberData(1).dVal = dVal;
    
    dataFileFiber = [dVal ' Fiber Data'];
    save([dataPath dataFileFiber],'fiberData'); 
end