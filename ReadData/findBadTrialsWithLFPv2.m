% We check, for each trial, whether there is a value greater than 
% 1. threshold times the std dev. 
% 2. maxLimit
% If yes, that trial is marked as a bad trial

function [allBadTrials,badTrials] = findBadTrialsWithLFPv2(monkeyName,expDate,protocolName,folderSourceString,gridType,checkTheseElectrodes,processAllElectrodes,threshold,maxLimit,showElectrodes,minLimit,saveDataFlag,checkPeriod,rejectTolerance)

if ~exist('checkTheseElectrodes','var');     checkTheseElectrodes = [33 12 80 63 44];   end
if ~exist('processAllElectrodes','var');     processAllElectrodes = 0;                  end
if ~exist('folderSourceString','var');       folderSourceString = 'G:';                 end
if ~exist('threshold','var');                threshold = 6;                             end
if ~exist('minLimit','var');                 minLimit = -2000;                          end
if ~exist('saveDataFlag','var');             saveDataFlag = 1;                          end
if ~exist('checkPeriod','var');             checkPeriod = [-0.7 0.8];                   end

folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName);
folderSegment = fullfile(folderName,'segmentedData');

load(fullfile(folderSegment,'LFP','lfpInfo.mat'));

if processAllElectrodes % compute bad trials for all the saved electrodes
    numElectrodes = length(analogChannelsStored);
else % compute bad trials for only the electrodes mentioned
    numElectrodes = length(checkTheseElectrodes);
end

allBadTrials = cell(1,numElectrodes);
nameElec = cell(1,numElectrodes);

for i=1:numElectrodes
    electrodeNum=analogChannelsStored(i);
    load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeNum) '.mat']));
    
    disp(['Processing electrode: ' num2str(electrodeNum)]);
    nameElec{i} = ['elec' num2str(electrodeNum)];
    disp(nameElec{i});
    
    analogDataSegment = analogData;
    % determine indices corresponding to the check period
    checkPeriodIndices = timeVals>=checkPeriod(1) & timeVals<=checkPeriod(2);
    
    analogData = analogData(:,checkPeriodIndices);
    % subtract dc
    analogData = analogData - repmat(mean(analogData,2),1,size(analogData,2));
    
    numTrials = size(analogData,1); %#ok<*NODEF>
    meanData = mean(analogData,2)';
    stdData  = std(analogData,[],2)';
    maxData  = max(analogData,[],2)';
    minData  = min(analogData,[],2)';
    
    clear tmpBadTrials tmpBadTrials2
    tmpBadTrials = unique([find(maxData > meanData + threshold * stdData) find(minData < meanData - threshold * stdData)]);
    tmpBadTrials2 = unique(find(maxData > maxLimit));
    tmpBadTrials3 = unique(find(minData < minLimit));
    allBadTrials{i} = unique([tmpBadTrials tmpBadTrials2 tmpBadTrials3]);
end

numElectrodes = length(checkTheseElectrodes); % check the list for these electrodes only to generate the overall badTrials list
j = find(analogChannelsStored==checkTheseElectrodes(1));
badTrials=allBadTrials{j};
for i=1:numElectrodes
    j = find(analogChannelsStored==checkTheseElectrodes(i));
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
            if ~isempty(find(allBadTrials{i} == n, 1))
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
    j = find(analogChannelsStored==checkTheseElectrodes(i));
    if length(allBadTrials{i}) ~= length(badTrials)
        disp(['Bad trials for electrode ' num2str(checkTheseElectrodes(i)) ': ' num2str(length(allBadTrials{j}))]);
    else
        disp(['Bad trials for electrode ' num2str(checkTheseElectrodes(i)) ': common bad trials only (' num2str(length(badTrials)) ')']);
    end
end

if saveDataFlag
    disp(['Saving ' num2str(length(badTrials)) ' bad trials']);
    save(fullfile(folderSegment,'badTrials.mat'),'badTrials','checkTheseElectrodes','threshold','maxLimit','checkPeriod','allBadTrials','nameElec');
else
    disp('Bad trials will not be saved..');
end

load(fullfile(folderSegment,'LFP','lfpInfo.mat'));

lengthShowElectrodes = length(showElectrodes);
if ~isempty(showElectrodes)
    for i=1:lengthShowElectrodes
        figure;
        if lengthShowElectrodes>1
            subplot(lengthShowElectrodes,1,i);
        else
            subplot(2,1,1);
        end
        channelNum = showElectrodes(i);

        clear signal analogData
        load(fullfile(folderSegment,'LFP',['elec' num2str(channelNum) '.mat']));
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
        
        if lengthShowElectrodes==1
            subplot(2,1,2);
            plot(timeVals,analogDataSegment(setdiff(1:numTrials,badTrials),:),'color','k');
            hold on;
            j = analogChannelsStored==channelNum;
            if ~isempty(allBadTrials{j})
                plot(timeVals,analogDataSegment(allBadTrials{j},:),'color','r');
            end
            axis tight;
        end
    end
end