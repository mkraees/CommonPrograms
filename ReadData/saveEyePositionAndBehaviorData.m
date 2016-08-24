% The old saveLLData program contained the following functions from the old format
% 1. saveLLDataSRC and saveLLDataGRF
% 2. getStimResultsLLSRC and getStimResultsLLGRF
% 3. saveEyeAndBehaviorDataSRC and saveEyeAndBehaviorDataGRF
% 4. getEyePositionAndBehavioralDataSRC and getEyePositionAndBehavioralDataGRF

% 5. saveEyeDataInDegSRC, saveEyeDataInDeg (GRF suffix added),
% saveEyeDataStimPosGRF and getEyeDataStimPosGRF

% Now saveLLData only contains the old % 1. saveLLDataSRC and
% saveLLDataGRF, while saveEyePositionAndBehaviorData contains the
% remaining four (2-5 above)

% Removed eyeRangeMS, eyeRangeLongMS and maxStimPos. Those are read from
% the LL file directly.

function saveEyePositionAndBehaviorData(subjectName,expDate,protocolName,folderSourceString,gridType,FsEye)

if ~exist('FsEye','var');        FsEye=200;                             end

folderName    = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');
makeDirectory(folderExtract);

if strncmpi(protocolName,'SRC',3) % SRC
    
    [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataSRC(subjectName,expDate,protocolName,folderSourceString,FsEye); 
    save(fullfile(folderExtract,'BehaviorData.mat'),'allTrials','goodTrials','stimData');
    save(fullfile(folderExtract,'EyeData.mat'),'eyeData','eyeRangeMS');
    
    saveEyeDataInDegSRC(subjectName,expDate,protocolName,folderSourceString,gridType);
else

    if strncmpi(protocolName,'GRF',3)
        [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataGRF(subjectName,expDate,protocolName,folderSourceString,FsEye); %#ok<*ASGLU,*NASGU>
    elseif strncmpi(protocolName,'CRS',3)
        [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataCRS(subjectName,expDate,protocolName,folderSourceString,FsEye); %#ok<*ASGLU,*NASGU>
    end
    save(fullfile(folderExtract,'BehaviorData.mat'),'allTrials','goodTrials','stimData');
    save(fullfile(folderExtract,'EyeData.mat'),'eyeData','eyeRangeMS');
    
   saveEyeDataInDegGRF(subjectName,expDate,protocolName,folderSourceString,gridType,FsEye);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataSRC(subjectName,expDate,protocolName,folderSourceString,Fs) %#ok<*DEFNU>

if ~exist('Fs','var');                   Fs = 200;                      end    % Eye position sampled at 200 Hz.

datFileName = fullfile(folderSourceString,'data','rawData',[subjectName expDate],[subjectName expDate protocolName '.dat']);

% Get Lablib data
header = readLLFile('i',datFileName);

if isfield(header,'cueDurationMS')
    cueDurationMS    = header.cueDurationMS.data;
    precueDurationMS = header.precueDurationMS.data;
    precueJitterPC   = header.precueJitterPC.data;
else
    cueDurationMS = 0;
    precueDurationMS = 500;
    precueJitterPC   = 0;
end
minFixationDurationMS = cueDurationMS + round((1-precueJitterPC/100) * precueDurationMS);

minStimDurationMS = round((1-header.stimJitterPC.data/100) * header.stimDurationMS.data);
minInterstimDurationMS = round((1-header.interstimJitterPC.data/100) * header.interstimMS.data);

eyeRangeMS = [-min(minFixationDurationMS,minInterstimDurationMS)+1000/Fs minStimDurationMS-1000/Fs]; % Around each stimulus onset, data should be available for this range. 2 samples are reduced on each range because sometimes there are minor timing shifts and 2 samples may not be collected.
eyeRangePos = eyeRangeMS*Fs/1000;

% Stimulus properties
numTrials = header.numberOfTrials;
stimNumber=1;
correctIndex=1;
trialEndIndex=1;

for i=1:numTrials
    %disp(i);
    clear trials
    trials = readLLFile('t',i);
    
    if isfield(trials,'trialEnd')
        allTrials.trialEnded(i) = 1;
        allTrials.catchTrials(trialEndIndex) = trials.trial.data.catchTrial;
        allTrials.instructTrials(trialEndIndex) = trials.trial.data.instructTrial;
        allTrials.trialCertify(trialEndIndex) = trials.trialCertify.data;
        allTrials.targetPosAllTrials(trialEndIndex) = trials.trial.data.targetIndex+1;
        allTrials.eotCodes(trialEndIndex) = trials.trialEnd.data;
        
        allTrials.fixWindowSize(trialEndIndex) = trials.fixWindowData.data.windowDeg.size.width;
        %allTrials.respWindowSize(trialEndIndex) = trials.respWindowData.data.windowDeg.size.width; % instead of responseWindowData
        allTrials.certifiedNonInstruction(trialEndIndex) = (allTrials.instructTrials(trialEndIndex)==0)*(allTrials.trialCertify(trialEndIndex)==0);

        if (allTrials.eotCodes(trialEndIndex)==0) &&  (allTrials.certifiedNonInstruction(trialEndIndex)==1) ...
                && (allTrials.catchTrials(trialEndIndex)==0) % Work on only Correct Trials, which are not instruction, catch or uncertified trials
            
            % Get Eye Data
            eyeX = trials.eyeXData.data;
            eyeY = trials.eyeYData.data;
            % eyeStartTime = trials.eyeXData.timeMS(1);  % This is wrong.
            % The eye data is synchronized with trialStartTime.
            eyeStartTime = trials.trialStart.timeMS;
            eyeAllTimes = eyeStartTime + (0:(length(eyeX)-1))*(1000/Fs);
            
            stimOnTimes  = [trials.stimulusOn.timeMS];
            numStimuli = allTrials.targetPosAllTrials(trialEndIndex); %=length(stimOnTimes)/3;
            
            goodTrials.targetPos(correctIndex) = numStimuli;
            goodTrials.targetTime(correctIndex) = stimOnTimes(end);
            if isfield(trials,'fixate') % [Vinay] - This variable/field is absent if task does not require fixation
                goodTrials.fixateMS(correctIndex) = trials.fixate.timeMS;
            end
            goodTrials.fixonMS(correctIndex) = trials.fixOn.timeMS;
            goodTrials.stimOnTimes{correctIndex} = stimOnTimes;
            
            clear stimType
            stimType = ([trials.stimDesc.data.type0]) .* ([trials.stimDesc.data.type1]);
            goodTrials.stimType{correctIndex} = stimType;
            
            for j=1:numStimuli

                if (stimType(j)==9) || (stimType(j)==1)  % Frontpad or Valid
                    
                    stimTime = stimOnTimes(j);
                    stp=find(eyeAllTimes>=stimTime, 1 );
                    
%                     stimData.stimOnsetTimeFromFixate(stimNumber) = stimTime-trials.fixate.timeMS;
                    if isfield(trials,'fixate') % [Vinay] - This variable/field is absent if task does not require fixation
                        stimData.stimOnsetTimeFromFixate(stimNumber) = stimTime-trials.fixate.timeMS;
                    end
                    stimData.stimPos(stimNumber) = j;
                    stimData.stimType(stimNumber) = stimType(j);
                    
                    if (stimType(j)==9) % First stimulus may not have sufficient baseline
                        eyeData(stimNumber).eyePosDataX = eyeX(stp:stp+eyeRangePos(2)-1); %#ok<*AGROW>
                        eyeData(stimNumber).eyePosDataY = eyeY(stp:stp+eyeRangePos(2)-1);
                    else
                        eyeData(stimNumber).eyePosDataX = eyeX(stp+eyeRangePos(1):stp+eyeRangePos(2)-1);
                        eyeData(stimNumber).eyePosDataY = eyeY(stp+eyeRangePos(1):stp+eyeRangePos(2)-1);  
                    end
                    
                    eyeData(stimNumber).eyeCal = trials.eyeCalibrationData.data.cal;
                    stimNumber=stimNumber+1;
                end
            end
            
            correctIndex=correctIndex+1;
        end
        trialEndIndex=trialEndIndex+1;
    end
end
end
function [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataGRF(subjectName,expDate,protocolName,folderSourceString,Fs)

if ~exist('Fs','var');                  Fs = 200;                       end    % Eye position sampled at 200 Hz.

datFileName = fullfile(folderSourceString,'data','rawData',[subjectName expDate],[subjectName expDate protocolName '.dat']);

% Get Lablib data
header = readLLFile('i',datFileName);

minFixationDurationMS = round((1-header.behaviorSetting.data.fixateJitterPC/100) * header.behaviorSetting.data.fixateMS);
stimDurationMS = header.mapStimDurationMS.data;
interStimDurationMS = header.mapInterstimDurationMS.data;

eyeRangeMS = [-min(minFixationDurationMS,interStimDurationMS)+1000/Fs stimDurationMS-1000/Fs]; % Around each stimulus onset, data should be available for this range. 2 samples are reduced on each range because sometimes there are minor timing shifts and 2 samples may not be collected.
eyeRangePos = eyeRangeMS*Fs/1000;

% Stimulus properties
numTrials = header.numberOfTrials;
stimNumber=1;
correctIndex=1;
trialEndIndex=1;

for i=1:numTrials
    disp(['Behavior: trial ' num2str(i) ' of ' num2str(numTrials)]);
    clear trials
    trials = readLLFile('t',i);
    
    if isfield(trials,'trialEnd')
        allTrials.trialEnded(i) = 1;
        allTrials.catchTrials(trialEndIndex) = trials.trial.data.catchTrial;
        allTrials.instructTrials(trialEndIndex) = trials.trial.data.instructTrial;
        allTrials.trialCertify(trialEndIndex) = trials.trialCertify.data;
        allTrials.targetPosAllTrials(trialEndIndex) = trials.trial.data.targetIndex+1;
        allTrials.eotCodes(trialEndIndex) = trials.trialEnd.data;
        
        allTrials.fixWindowSize(trialEndIndex) = trials.fixWindowData.data.windowDeg.size.width;
        allTrials.respWindowSize(trialEndIndex) = trials.responseWindowData.data.windowDeg.size.width;
        allTrials.certifiedNonInstruction(trialEndIndex) = (allTrials.instructTrials(trialEndIndex)==0)*(allTrials.trialCertify(trialEndIndex)==0);

        if (allTrials.eotCodes(trialEndIndex)==0) &&  (allTrials.certifiedNonInstruction(trialEndIndex)==1)
                %&& (allTrials.catchTrials(trialEndIndex)==0) % Work on only Correct Trials, which are not instruction or uncertified trials. Include catch trials
            
            isCatchTrial = (allTrials.catchTrials(trialEndIndex)==1);
            
            % Get Eye Data
            if isfield(trials,'eyeXData')
                eyeX = trials.eyeXData.data;
                eyeY = trials.eyeYData.data;
            elseif isfield(trials,'eyeRXData')
                eyeX = trials.eyeRXData.data;
                eyeY = trials.eyeRYData.data;
            elseif isfield(trials,'eyeLXData')
                eyeX = trials.eyeLXData.data;
                eyeY = trials.eyeLYData.data;
            end
            
            % eyeStartTime = trials.eyeXData.timeMS(1);  % This is wrong.
            % The eye data is synchronized with trialStartTime.
            % eyeStartTime = trials.trialStart.timeMS;
            
            % Not any more. Now after trialStart, we sleep for sometime to
            % send long digital pulses. Now we use the start of eye
            % calibration as the onset time.
            eyeStartTime = trials.eyeLeftCalibrationData.timeMS;
            eyeAllTimes = eyeStartTime + (0:(length(eyeX)-1))*(1000/Fs);
            
            stimOnTimes  = [trials.stimulusOnTime.timeMS];
            numStimuli = length(stimOnTimes)/3; % = allTrials.targetPosAllTrials(trialEndIndex); %=length(stimOnTimes)/3;
            
            goodTrials.targetPos(correctIndex) = numStimuli;
            goodTrials.targetTime(correctIndex) = stimOnTimes(end);
            if isfield(trials,'fixate') % [Vinay] - This variable/field is absent if task does not require fixation
                goodTrials.fixateMS(correctIndex) = trials.fixate.timeMS;
            end
            goodTrials.fixonMS(correctIndex) = trials.fixOn.timeMS;
            goodTrials.stimOnTimes{correctIndex} = stimOnTimes;
            
            % Find position of Gabor1
            gaborPos = find([trials.stimDesc.data.gaborIndex]==1); % could be 4 gabors for GRF protocol
            
%             if isCatchTrial
%                 stimEndIndex = numStimuli;    % Take the last stimulus because it is still a valid stimulus
%             else
%                 stimEndIndex = numStimuli-1;  % Don't take the last one because it is the target
%             end
            
            if trials.stimDesc.data(1).stimType == 0 % If the target gabor is null then don't check for catch trial and consider all stimuli as valid for analyses
                stimEndIndex = numStimuli;
            else
                if isCatchTrial
                    stimEndIndex = numStimuli;    % Take the last stimulus because it is still a valid stimulus
                else
                    stimEndIndex = numStimuli-1;  % Don't take the last one because it is the target
                end
            end
                
            if stimEndIndex>0  % At least one stimulus
                for j=1:stimEndIndex
                    stimTime = stimOnTimes(gaborPos(j));
                    stp=find(eyeAllTimes>=stimTime,1);
                    
                    if isfield(trials,'fixate') % [Vinay] - This variable/field is absent if task does not require fixation
                        stimData.stimOnsetTimeFromFixate(stimNumber) = stimTime-trials.fixate.timeMS;
                    end
                    stimData.stimPos(stimNumber) = j;
                    
                    startingPos = max(1,stp+eyeRangePos(1));
                    endingPos   = min(stp+eyeRangePos(2)-1,length(eyeX));
     
                    eyeData(stimNumber).eyePosDataX = eyeX(startingPos:endingPos);
                    eyeData(stimNumber).eyePosDataY = eyeY(startingPos:endingPos);
                    
                    if isfield(trials,'eyeXData')
                        eyeData(stimNumber).eyeCal = trials.eyeCalibrationData.data.cal;
                    elseif isfield(trials,'eyeRXData')
                        eyeData(stimNumber).eyeCal = trials.eyeRightCalibrationData.data.cal;
                    elseif isfield(trials,'eyeLXData')
                        eyeData(stimNumber).eyeCal = trials.eyeLeftCalibrationData.data.cal;
                    end
                    
                    stimNumber=stimNumber+1;
                end
            end
            correctIndex=correctIndex+1;
        end
        trialEndIndex=trialEndIndex+1;
    end
end
end
% Additional analysis
function saveEyeDataInDegSRC(subjectName,expDate,protocolName,folderSourceString,gridType)
% The difference between saveEyeDataInDeg and saveEyeDataInDegSRC is that
% the frontPad stimuli are also saved in SRC.

folderName    = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');

clear eyeData 
load(fullfile(folderExtract,'EyeData.mat'));

[eyeDataDegX,eyeDataDegY] = convertEyeDataToDeg(eyeData,1);

for i=1:length(eyeDataDegX)
    lengthEyeSignal = size(eyeDataDegX{i},1);
    eyeSpeedX{i} = [eyeDataDegX{i}(2:lengthEyeSignal)-eyeDataDegX{i}(1:lengthEyeSignal-1);0];
    eyeSpeedY{i} = [eyeDataDegY{i}(2:lengthEyeSignal)-eyeDataDegY{i}(1:lengthEyeSignal-1);0];
end

folderSave = fullfile(folderName,'segmentedData','eyeData');
makeDirectory(folderSave);
save(fullfile(folderSave,'eyeDataDeg.mat'),'eyeDataDegX','eyeDataDegY');
save(fullfile(folderSave,'eyeSpeed.mat'),'eyeSpeedX','eyeSpeedY');

end
function saveEyeDataInDegGRF(subjectName,expDate,protocolName,folderSourceString,gridType,FsEye)

folderName    = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');

clear eyeData 
load(fullfile(folderExtract,'EyeData.mat'));

clear goodStimNums
load(fullfile(folderExtract,'goodStimNums.mat'));

if exist(fullfile(folderExtract,'validStimAfterTarget.mat'),'file')
    load(fullfile(folderExtract,'validStimAfterTarget.mat'));
    disp(['Removing ' num2str(length(validStimuliAfterTarget)) ' stimuli from goodStimNums']);
    goodStimNums(validStimuliAfterTarget)=[];
end

clear stimResults
load(fullfile(folderExtract,'stimResults.mat'));

goodStimPos = stimResults.stimPosition(goodStimNums);

% all stimPostions greater than 1
% useTheseStims = find(goodStimPos>1);
useTheseStims = find(goodStimPos>0); % Use all stimPositions, including 1

% [Vinay] - Decide whether to return a cell array for eyeDataDegX/Y
% depending on eyePosDataX/Y's length across trials. If the
% lengths are not same then construct eyeDataDegX/Y as cell arrays
% Non-cell arrays, in this case, would flag dimension mismatch
for i = 1:length(eyeData)
    lenEyePos(i) = length(eyeData(i).eyePosDataX);
end

returnCellArray = 0;
if ~isequal(diff(lenEyePos),zeros(size(lenEyePos)))
    returnCellArray = 1;
end
%---

[eyeDataDegX,eyeDataDegY] = convertEyeDataToDeg(eyeData(useTheseStims),returnCellArray); % [Vinay] - pass returnCellArray (instead of default 0)
folderSave = fullfile(folderName,'segmentedData','eyeData');
makeDirectory(folderSave);
save(fullfile(folderSave,'eyeDataDeg.mat'),'eyeDataDegX','eyeDataDegY');

% lengthEyeSignal = size(eyeDataDegX,2);
% eyeSpeedX = [eyeDataDegX(:,2:lengthEyeSignal)-eyeDataDegX(:,1:lengthEyeSignal-1) zeros(size(eyeDataDegX,1),1)];
% eyeSpeedY = [eyeDataDegY(:,2:lengthEyeSignal)-eyeDataDegY(:,1:lengthEyeSignal-1) zeros(size(eyeDataDegY,1),1)];
% save([folderSave 'eyeSpeed.mat'],'eyeSpeedX','eyeSpeedY');

% More data saved for GRF protocol
[eyeXAllPos,eyeYAllPos,xs,durationsMS] = getEyeDataStimPosGRF(subjectName,expDate,protocolName,folderSourceString,FsEye);
save(fullfile(folderSave,'EyeDataStimPos.mat'),'eyeXAllPos','eyeYAllPos','xs','durationsMS');
end
function [eyeXAllPos,eyeYAllPos,xs,durationsMS] = getEyeDataStimPosGRF(subjectName,expDate,protocolName,folderSourceString,FsEye)

intervalTimeMS=1000/FsEye;
datFileName = fullfile(folderSourceString,'data','rawData',[subjectName expDate],[subjectName expDate protocolName '.dat']);

% Get Lablib data
header = readLLFile('i',datFileName);

minFixationDurationMS = round((1-header.behaviorSetting.data.fixateJitterPC/100) * header.behaviorSetting.data.fixateMS);
stimDurationMS = header.mapStimDurationMS.data;
interStimDurationMS = header.mapInterstimDurationMS.data;
% <<<<<<< HEAD
% % maxStimPos = ceil(header.maxTargetTimeMS.data)/(stimDurationMS+interStimDurationMS) +1;
% maxStimPos = ceil((header.maxTargetTimeMS.data)/(stimDurationMS+interStimDurationMS)) +1; % [Vinay] - ceil should be for the entire expression?
% =======
responseTimeMS = header.responseTimeMS.data;
maxStimPos = ceil(header.maxTargetTimeMS.data + responseTimeMS)/(stimDurationMS+interStimDurationMS) +1;
% >>>>>>> master

durationsMS.minFixationDurationMS = minFixationDurationMS;
durationsMS.interStimDurationMS = interStimDurationMS;
durationsMS.stimDurationMS = stimDurationMS;

for j=1:maxStimPos                                  
    eyeXAllPos{j}=[];
    eyeYAllPos{j}=[];
    xs{j}=[];
end
 
numTrials = header.numberOfTrials;

for i=1:numTrials
    trial = readLLFile('t',i);
    disp(['Extended EyeData: trial ' num2str(i) ' of ' num2str(numTrials)]);
    % Work on only Correct Trials, which are not instruction or uncertified
    % trials. Include catch trials
    if (trial.trialEnd.data == 0) && (trial.trial.data.instructTrial==0) && ...
            (trial.trialCertify.data==0) %&& (trial.trial.data.catchTrial==0)
        
        isCatchTrial = (trial.trial.data.catchTrial==1);
        
        % get eye data
        clear eX eY cal
        if isfield(trial,'eyeXData')
            eX = trial.eyeXData.data';
            eY = trial.eyeYData.data';
            cal=trial.eyeCalibrationData.data.cal;
            timeStartMS = trial.eyeCalibrationData.timeMS;
        elseif isfield(trial,'eyeRXData')
            eX = trial.eyeRXData.data';
            eY = trial.eyeRYData.data';
            cal=trial.eyeRightCalibrationData.data.cal;
            timeStartMS = trial.eyeRightCalibrationData.timeMS;
        elseif isfield(trial,'eyeLXData')
            eX = trial.eyeLXData.data';
            eY = trial.eyeLYData.data';
            cal=trial.eyeLeftCalibrationData.data.cal;
            timeStartMS = trial.eyeLeftCalibrationData.timeMS;
        end
        
        if isCatchTrial
            numUsefulStim = trial.trial.data.numStim; % these are the useful stimuli, including target.
        else
            numUsefulStim = trial.trial.data.targetIndex; % these are the useful stimuli, excluding target.
        end
        
        stimOnTimes = trial.stimulusOnTime.timeMS;
        gaborPos = find([trial.stimDesc.data.gaborIndex]==1);
        
        if numUsefulStim>0
            for j=1:numUsefulStim
                stimOnsetPos = ceil((stimOnTimes(gaborPos(j)) - timeStartMS)/intervalTimeMS);
                
                stp = -(minFixationDurationMS + (j-1)*(stimDurationMS+interStimDurationMS))/intervalTimeMS + 1;                  
                edp = stimDurationMS/intervalTimeMS - 1;
                list = stp:edp;
                
                % [Vinay] - shift list if there are not enough data points
                % in eX/Y
                wrtStimOnsetPos = stimOnsetPos+list;
                shiftList = 0;
                if(wrtStimOnsetPos(end)-length(eX)>0)
                    shiftList = wrtStimOnsetPos(end)-length(eX);
                end
                
%                 eXshort = eX(stimOnsetPos+list);
%                 eYshort = eY(stimOnsetPos+list);
                
%                 eXshort = eX(stimOnsetPos+list-shiftList);
%                 eYshort = eY(stimOnsetPos+list-shiftList);
                if(wrtStimOnsetPos(1)-length(eX)>0) % No valid eye data available (Why does this happen?)
                    listLen = length(list);
                    eXshort = [];
                    eYshort = [];
                    eXshortDeg = [];
                    eYshortDeg = [];
                    eXshort(1:listLen) = 0;
                    eYshort(1:listLen) = 0;
                    eXshortDeg(1:listLen) = 0;
                    eYshortDeg(1:listLen) = 0;
                    disp(['No eye data! trial no. :' num2str(i)]);
                else
                    listLen = length(list);
                    list(end:-1:end-shiftList+1) = []; % to read only till the end of eX or eY and not exceed their indices
                    eXshort = eX(stimOnsetPos+list);
                    eYshort = eY(stimOnsetPos+list);
                    % append zeros to these if they are not of length=listLen,
                    % else there is index mismatch
                    eXshort(end+1:listLen) = 0;
                    eYshort(end+1:listLen) = 0;

                    eXshortDeg = cal.m11*eXshort + cal.m21 * eYshort + cal.tX;
                    eYshortDeg = cal.m12*eXshort + cal.m22 * eYshort + cal.tY;
                end
                
                eyeXAllPos{j} = cat(1,eyeXAllPos{j},eXshortDeg);
                eyeYAllPos{j} = cat(1,eyeYAllPos{j},eYshortDeg);
                
                if isempty(xs{j})
                    xs{j} = list*intervalTimeMS + (j-1)*(stimDurationMS+interStimDurationMS);
                end
            end
        end
    end
end
end