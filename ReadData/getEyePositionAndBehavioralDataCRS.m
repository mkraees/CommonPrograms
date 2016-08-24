% Removed from saveEyePositionAndBehaviorData.m and stored separately for 
% CRS
% 24 August 2016
% *************************************************************************
% Read and store eye position and behaviour information from Lablib data 
% file
% originally modified from the function: 
% 'getEyePositionAndBehavioralDataGRF' in saveEyePositionAndBehaviorData.m
% 
%**************************************************************************

function [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataCRS(subjectName,expDate,protocolName,folderSourceString,Fs)

if ~exist('Fs','var');                  Fs = 200;                       end    % Eye position sampled at 200 Hz.

datFileName = fullfile(folderSourceString,'data','rawData',[subjectName expDate],[subjectName expDate protocolName '.dat']);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
datFileInfoName = fullfile(folderSourceString,'data','rawData',[subjectName expDate],[subjectName expDate protocolName '_fileInfo.mat']);
removeFile=0; % [Vinay] - added this loop to delete the already generated files (due to which MATLAB keeps throwing an error)
if removeFile
    if exist(datFileInfoName,'file')
        delete(datFileInfoName);
    end
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
            % [Vinay] - for CRS the number of gabors is 4, hence the
            % denominator changes to 4 from 3 in the next line
            numStimuli = length(stimOnTimes)/4; % = allTrials.targetPosAllTrials(trialEndIndex); %=length(stimOnTimes)/3;
                
            goodTrials.targetPos(correctIndex) = numStimuli;
            goodTrials.targetTime(correctIndex) = stimOnTimes(end);
            if isfield(trials,'fixate') % [Vinay] - This variable/field is absent if task does not require fixation
                goodTrials.fixateMS(correctIndex) = trials.fixate.timeMS;
            end
            goodTrials.fixonMS(correctIndex) = trials.fixOn.timeMS;
            goodTrials.stimOnTimes{correctIndex} = stimOnTimes;
            
            % Find position of Gabor1
            gaborPos = find([trials.stimDesc.data.gaborIndex]==1); % could be 4 gabors for GRF protocol
            
            if trials.stimDesc.data(1).stimType == 0 % [Vinay] - If the target gabor is null then don't check for catch trial and consider all stimuli as valid for analyses
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