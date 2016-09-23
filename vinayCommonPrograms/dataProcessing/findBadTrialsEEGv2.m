% find bad trials for each electrode for EEG
% modified from findBadTrialsWithLFP and runFindBadTrials
%
% Vinay Shirhatti, 12 Oct 2014
%
% 24 March 2015, modified from findBadTrialsEEG in
% GammaStimDiscontinuityProject
% 21 Sep 2016, modified as per findBadTrialsLFPv2
%==========================================================================


function [allBadTrials, badTrials, nameElec] = findBadTrialsEEGv2(subjectName,expDate,protocolName,folderSourceString,gridType,...
    checkTheseElectrodes,processAllElectrodes,threshold,maxLimit,minLimit,showElectrodes,saveDataFlag,checkPeriod,rejectTolerance,showTrials)

if ~exist('processAllElectrodes','var');    processAllElectrodes = 0;      end
if ~exist('folderSourceString','var')       folderSourceString = 'K:\';    end
if ~exist('threshold','var')                threshold = 6;                 end
if ~exist('minLimit','var')                 minLimit = -300;               end
if ~exist('maxLimit','var')                 maxLimit = 300;                end
if ~exist('saveDataFlag','var');            saveDataFlag = 1;              end
if ~exist('checkPeriod','var');             checkPeriod = [-0.7 0.8];      end
if ~exist('rejectTolerance','var');         rejectTolerance = 1;           end

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');

load(fullfile(folderLFP,'lfpInfo'));

if ~exist('checkTheseElectrodes','var')
    checkTheseElectrodes = analogChannelsStored;   
end

if processAllElectrodes % compute bad trials for all the saved electrodes
    numElectrodes = length(analogChannelsStored);
else % compute bad trials for only the electrodes mentioned
    numElectrodes = length(checkTheseElectrodes);
end

allBadTrials = cell(1,numElectrodes);
nameElec = cell(1,numElectrodes);

for i=1:numElectrodes
    clear analogData
    
    if processAllElectrodes
        electrodeNum=analogChannelsStored(i); % changed from checkTheseElectrodes
        % to calculate bad trials for each electrode irrespective of the
        % electrodes to be checked
    else
        electrodeNum=checkTheseElectrodes(i);
    end
    load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeNum) '.mat']));
    
    disp(['Processing electrode: ' num2str(electrodeNum)]);
    nameElec{i} = ['elec' num2str(electrodeNum)];
    disp(nameElec{i});
    
    % Set the limits higher for frontal electrodes 1,2 & 61(for eyeblinks)
    if (electrodeNum == 1 || electrodeNum == 2 || electrodeNum == 61)
        maxLimitElec = 300; minLimitElec = -300;
    else
        maxLimitElec = maxLimit; minLimitElec = minLimit;
    end
    
    analogDataSegment = analogData;
    % determine indices corresponding to the check period
    checkPeriodIndices = timeVals>=checkPeriod(1) & timeVals<=checkPeriod(2);
    
    analogData = analogData(:,checkPeriodIndices);
    % subtract dc
    analogData = analogData - repmat(mean(analogData,2),1,size(analogData,2));
    
    if showTrials
        hAllTrials = figure(11);
        subplot(8,8,i); plot(timeVals,analogData);title(['elec' num2str(i)]);
        axis('tight');
    end
    
    numTrials = size(analogData,1); %#ok<*NODEF>
    meanData = mean(analogData,2)';
    stdData  = std(analogData,[],2)';
    maxData  = max(analogData,[],2)';
    minData  = min(analogData,[],2)';
    
    clear tmpBadTrials1 tmpBadTrials2 tmpBadTrials3 tmpBadTrials4 
    tmpBadTrials1 = unique([find(maxData > meanData + threshold * stdData) find(minData < meanData - threshold * stdData)]);
    % Vinay - Ideally set maxLimit and minLimit for the below criteria to
    % be quite high/low so that only the extremely noisy trials are
    % rejected by these
    tmpBadTrials2 = unique(find(maxData > maxLimitElec));
    tmpBadTrials3 = unique(find(minData < minLimitElec));
    
    % Vinay - Set another criterion based on the deviation of each trial
    % from the mean trial signal. This way even if the actual signal shows
    % high deflections, if those deflections are consistent then the trials
    % are not rejected purely on the max/min criteria above
    meanTrialData = mean(analogData,1); % mean trial trace
    stdTrialData = std(analogData,[],1); % std across trials
%     maxTrialData = max(analogData,[],1);
%     minTrialData = min(analogData,[],1);
    tDplus = (meanTrialData + (threshold)*stdTrialData); % upper boundary/criterion
    tDminus = (meanTrialData - (threshold)*stdTrialData); % lower boundary/criterion
    
    % Check for trials that cross these boundaries and mark them as bad
    tBoolTrials = zeros(1,numTrials);
    for tr = 1:numTrials
        
        trialDeviationHigh = tDplus - analogData(tr,:); % deviation from the upper boundary
        trialDeviationLow = analogData(tr,:) - tDminus; % deviation from the lower boundary
        
        tBool = zeros(size(meanTrialData,1),size(meanTrialData,2));
        tBool(trialDeviationHigh<0) = 1; % set if the upper boundary is crossed anywhere
        tBool1 = sum(tBool); % number of times the upper boundary was crossed
        
        tBool = zeros(size(meanTrialData,1),size(meanTrialData,2));
        tBool(trialDeviationLow<0) = 1; % set if the lower boundary is crossed anywhere
        tBool2 = sum(tBool); % number of times the lower boundary was crossed

        tBoolTrials(tr) = tBool1 + tBool2; % number of times the boundaries were crossed
    end
    
    tmpBadTrials4 = find(tBoolTrials>0);
    
    allBadTrials{electrodeNum} = unique([tmpBadTrials1 tmpBadTrials2 tmpBadTrials3 tmpBadTrials4]);
    
    goodTrials = 1:numTrials;
    goodTrials = setdiff(goodTrials,allBadTrials{electrodeNum});
    
    clear meanData maxData minData stdData meanTrialData stdTrialData maxLimitElec minLimitElec
    clear trialDeviationHigh trialDeviationLow tBool tBool1 tBool2 tBoolTrials tDminus tDplus
    
    if showTrials && ~isempty(goodTrials)
        hGoodTrials = figure(12);
        subplot(8,8,i); plot(timeVals,analogData(goodTrials,:));title(['elec' num2str(i)]);
        axis('tight');
    end
end

numElectrodes = length(checkTheseElectrodes); % check the list for these electrodes only to generate the overall badTrials list
j = checkTheseElectrodes(1);
badTrials=allBadTrials{j};
for i=1:numElectrodes
    j = checkTheseElectrodes(i);
    badTrials=intersect(badTrials,allBadTrials{j}); % in the previous case we took the union
end

disp(['total Trials: ' num2str(numTrials) ', bad trials: ' num2str(badTrials)]);

% [Vinay] - decide as per a tolerance for the percent of electrodes showing
% a particular stimulus as bad
if exist('rejectTolerance','var')
    badTrials = [];
    for n=1:numTrials
        trialCount=0;
        for i=1:numElectrodes
            j = checkTheseElectrodes(i);
            if ~isempty(find(allBadTrials{j} == n, 1))
                trialCount=trialCount+1;
            end
        end
        trialPercent = trialCount/numElectrodes;
        if trialPercent>=rejectTolerance
            badTrials = cat(1,badTrials,n);
        end
    end
end
%-----

for i=1:numElectrodes
    j = checkTheseElectrodes(i);
    if length(allBadTrials{j}) ~= length(badTrials)
        disp(['Bad trials for electrode ' num2str(checkTheseElectrodes(i)) ': ' num2str(length(allBadTrials{j}))]);
    else
        disp(['Bad trials for electrode ' num2str(checkTheseElectrodes(i)) ': common bad trials only (' num2str(length(badTrials)) ')']);
    end
end

if saveDataFlag
    disp(['Saving ' num2str(length(badTrials)) ' bad trials']);
    save(fullfile(folderSegment,'badTrials.mat'),'badTrials','checkTheseElectrodes','threshold','maxLimit','minLimit','checkPeriod','allBadTrials','nameElec','rejectTolerance');
else
    disp('Bad trials will not be saved..');
end

if ~isempty(showElectrodes)
    for i=1:length(showElectrodes)
        figure;
        subplot(2,1,1);
        channelNum = showElectrodes(i);

        clear signal analogData analogDataSegment
        load(fullfile(folderSegment,'LFP',['elec' num2str(channelNum) '.mat']));
        analogDataSegment = analogData;
        if numTrials<4000
            plot(timeVals,analogDataSegment(setdiff(1:numTrials,badTrials),:),'color','k');
            hold on;
        else
            disp('More than 4000 trials...');
        end
        if ~isempty(badTrials)
            plot(timeVals,analogDataSegment(badTrials,:),'color','g');
        end
        title(['electrode ' num2str(channelNum)]);
        axis tight;

        subplot(2,1,2);
        plot(timeVals,analogDataSegment(setdiff(1:numTrials,badTrials),:),'color','k');
        hold on;
        j = channelNum;
        if ~isempty(allBadTrials{j})
            plot(timeVals,analogDataSegment(allBadTrials{j},:),'color','r');
        end
        axis tight;
    end
end

end
