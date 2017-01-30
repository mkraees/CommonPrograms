
protocolIndices = 48:49;
[expDates,protocolNames,stimTypes,subjectNames] = allProtocolsHumanEEG;

for i=1:length(protocolIndices)
    subjectName = subjectNames{protocolIndices(i)}; 
    gridType = 'EEG';
    expDate = expDates{protocolIndices(i)}; 
    protocolName = protocolNames{protocolIndices(i)};
    FsEye=500; % 500Hz for eyelink
    timeStartFromBaseLineList(1) = -0.55; deltaTList(1) = 1.024; % in seconds
    timeStartFromBaseLineList(2) = -1.148; deltaTList(2) = 2.048;
    timeStartFromBaseLineList(3) = -1.096; deltaTList(3) = 4.096; % [Vinay] - for a longer stim time
    timeStartFromBaseLineList(4) = -0.848; deltaTList(4) = 2.048; % [Vinay] - 1000ms stim, 500ms interstim
    timeStartFromBaseLineList(5) = -1.596; deltaTList(5) = 4.096; % [Vinay] - for a longer stim time; Stim ON and OFF 1500ms each
    % folderSourceString = 'K:\';
    folderSourceString = '/media/vinay/SRLHD02M';
    electrodesToStore = [];
    ignoreTargetStimFlag=1;
    frameRate=100;
    convertToImageFlag = 0;
    reducedDigitalCodes = 1;
    deviceName = 'BP';
    stimType = stimTypes{protocolIndices(i)};
    type = stimType;
    deltaT = deltaTList(type);
    timeStartFromBaseLine = timeStartFromBaseLineList(type);
    folderExtract = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData');
    [digitalTimeStamps,digitalEvents]=extractDigitalDataBrainProducts(subjectName,expDate,protocolName,folderSourceString,gridType);
    saveDigitalData(digitalEvents,digitalTimeStamps,folderExtract);
    LLFileExistsFlag = saveLLData(subjectName,expDate,protocolName,folderSourceString,gridType);
    displayTSTEComparison(folderExtract);
    if strncmpi(protocolName,'GRF',3)
        [goodStimNums,goodStimTimes,activeSide]=extractDigitalDataGRFLL(folderExtract,ignoreTargetStimFlag,frameRate);
        getDisplayCombinationsGRF(folderExtract,goodStimNums);
    else
        [goodStimNums,goodStimTimes,activeSide]=extractDigitalDataCRSLL(folderExtract,ignoreTargetStimFlag,frameRate);
        getDisplayCombinationsCRS(folderExtract,goodStimNums);
    end
    saveEyePositionAndBehaviorData(subjectName,expDate,protocolName,folderSourceString,gridType,FsEye);
    getEEGDataBrainProducts(subjectName,expDate,protocolName,folderSourceString,gridType,goodStimTimes,timeStartFromBaseLine,deltaT);
end