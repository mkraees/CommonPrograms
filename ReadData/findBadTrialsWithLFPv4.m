% We check, for each trial, whether there is a value greater than 
% 1. threshold times the std dev. 
% 2. maxLimit
% If yes, that trial is marked as a bad trial
% 
% 07 July 2017, incorporating code for rejections based on statistics at
% every time point, as done for EEG
% modified from findBadTrialsWithLFPv3 and findBadTrialsEEGv2
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function [allBadTrials,badTrials] = findBadTrialsWithLFPv4(monkeyName,expDate,protocolName,folderSourceString,gridType,checkTheseElectrodes,processAllElectrodes,threshold,maxLimit,minLimit,showElectrodes,saveDataFlag,checkPeriod,rejectTolerance)

if ~exist('processAllElectrodes','var');     processAllElectrodes = 0;                  end
if ~exist('folderSourceString','var');       folderSourceString = 'G:';                 end
if ~exist('threshold','var');                threshold = 6;                             end
if ~exist('minLimit','var');                 minLimit = -2000;                          end
if ~exist('maxLimit','var');                 maxLimit = 1000;                           end
if ~exist('saveDataFlag','var');             saveDataFlag = 1;                          end
if ~exist('checkPeriod','var');              checkPeriod = [-0.7 0.8];                  end
if ~exist('rejectTolerance','var');          rejectTolerance = 1;                       end
if ~exist('showElectrodes','var');           showElectrodes = [];                       end

folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName);
folderSegment = fullfile(folderName,'segmentedData');

load(fullfile(folderSegment,'LFP','lfpInfo.mat'));

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
    maxLimitElec = maxLimit; minLimitElec = minLimit; % same limits for all electrodes unlike EEG where the threshold is higher for frontal electrodes
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
    
    clear meanData maxData minData stdData meanTrialData stdTrialData maxLimitElec minLimitElec
    clear trialDeviationHigh trialDeviationLow tBool tBool1 tBool2 tBoolTrials tDminus tDplus
end

numElectrodes = length(checkTheseElectrodes); % check the list for these electrodes only to generate the overall badTrials list
j = checkTheseElectrodes(1);
badTrials=allBadTrials{j};
for i=1:numElectrodes
    j = checkTheseElectrodes(i);
    badTrials=intersect(badTrials,allBadTrials{j}); % in the previous case we took the union
end
disp(['total Trials: ' num2str(numTrials) ', bad trials: ' num2str(badTrials)]);

%**************************************************************************
% [vinay] SOME ADDITIONAL CRITERIA for detecting bad repeats and electrodes
%--------------------------------------------------------------------------
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
            badTrials = cat(2,badTrials,n);
        end
    end
end
%-----

% [vinay] - statistically determine the repeats which are flagged as bad
% across a significant number of electrodes and reject them. Similarly tag
% those electrodes as bad which show a significantly higher number of
% badTrials as compared to the mean statistics across electrodes.

% look at the overall distribution of badTrials across electrodes
% and reject the repeats which are flagged as bad for >= (mean+std*threshold) 
% number of electrodes. (across trials)
% Similarly tag those electrodes as bad which show number of badTrials 
% >= (mean+std*threshold) (across electrodes)
allBadTrialsMatrix = zeros(length(allBadTrials),numTrials);
for i=1:length(allBadTrials)
    allBadTrialsMatrix(i,allBadTrials{i}) = 1;
end
marginalStim = sum(allBadTrialsMatrix(checkTheseElectrodes,:),1); % for each stimulus - no. of electrodes showing this stimulus as bad
marginalElectrodes = sum(allBadTrialsMatrix,2); % for each electrode - no. of stimuli flagged as bad

if processAllElectrodes && strcmpi(gridType,'Microelectrode')
    if strcmpi(monkeyName,'tutu') || strcmpi(monkeyName,'alpaH') || strcmpi(monkeyName,'kesariH')
        electrodesForMarginals = 1:81;
    else
        electrodesForMarginals = 1:96;
    end
else
    electrodesForMarginals = checkTheseElectrodes;
end

thresholdMarginal = threshold/2; % thresholdMarginal i.e. values> u+thresholdMarginal*sigma will be tagged as bad 
badElecsMarginalStats = [];
if mean(marginalElectrodes(electrodesForMarginals))>0
    badElecsMarginalStats = find(marginalElectrodes(electrodesForMarginals)>=(mean(marginalElectrodes(electrodesForMarginals))+std(marginalElectrodes(electrodesForMarginals))*thresholdMarginal));
end

% [vinay] 11 July 2017
% determine marginalstim after removing the badElecs and then determine
% badTrialsMarginalStats
marginalStim = sum(allBadTrialsMatrix(setdiff(checkTheseElectrodes,badElecsMarginalStats),:),1);
badTrialsMarginalStats = [];
if mean(marginalStim)>0
    badTrialsMarginalStats = find(marginalStim>=(mean(marginalStim)+std(marginalStim)*thresholdMarginal));
end

badTrials = unique(cat(2,badTrials,badTrialsMarginalStats));
badElecs = badElecsMarginalStats;

disp(['total Trials: ' num2str(numTrials) ', bad trials indices: ' num2str(badTrials)]);
disp(['total Elecs: ' num2str(length(marginalElectrodes)) ', bad elecs indices: ' num2str(badElecs')]);
%--------------------------------

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
    save(fullfile(folderSegment,'badTrials.mat'),'badTrials','checkTheseElectrodes','threshold','maxLimit','minLimit','checkPeriod','allBadTrials','nameElec','rejectTolerance','badElecs','badTrialsMarginalStats','allBadTrialsMatrix');
else
    disp('Bad trials will not be saved..');
end

load(fullfile(folderSegment,'LFP','lfpInfo.mat'));

lengthShowElectrodes = length(showElectrodes);
if ~isempty(showElectrodes)
    for i=1:lengthShowElectrodes
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

%**************************************************************************
% summary plot
%--------------------------------------------------------------------------
% allBadTrialsMatrix = zeros(length(allBadTrials),numTrials);
% for i=1:length(allBadTrials)
%     allBadTrialsMatrix(i,allBadTrials{i}) = 1;
% end

summaryFig = figure('name',[monkeyName expDate protocolName],'numbertitle','off');
h0 = subplot('position',[0.8 0.8 0.18 0.18]); set(h0,'visible','off');
text(0.05, 0.7, ['thresholds (uV): [' num2str(minLimit) ' ' num2str(maxLimit) ']'],'fontsize',11,'unit','normalized','parent',h0);
text(0.05, 0.55, ['checkPeriod (s): [' num2str(checkPeriod(1)) ' ' num2str(checkPeriod(2)) ']'],'fontsize',11,'unit','normalized','parent',h0);
text(0.05, 0.4, ['rejectTolerance : ' num2str(rejectTolerance)],'fontsize',11,'unit','normalized','parent',h0);
text(0.05, 0.25, ['threshold (*sig) : ' num2str(threshold)],'fontsize',11,'unit','normalized','parent',h0);
text(0.05, 0.1, ['thresholdMarginal : ' num2str(thresholdMarginal)],'fontsize',11,'unit','normalized','parent',h0);

h1 = getPlotHandles(1,1,[0.07 0.07 0.7 0.7]);
subplot(h1);
imagesc(1:numTrials,length(allBadTrials):-1:1,flipud(allBadTrialsMatrix),'parent',h1);
set(gca,'YDir','normal'); colormap(gray);
xlabel('# trial num','fontsize',15,'fontweight','bold');
ylabel('# electrode num','fontsize',15,'fontweight','bold');

h2 = getPlotHandles(1,1,[0.07 0.8 0.7 0.17]);
h3 = getPlotHandles(1,1,[0.8 0.07 0.18 0.7]);
subplot(h2); cla; set(h2,'nextplot','add');
stem(h2,1:numTrials,sum(allBadTrialsMatrix,1)); axis('tight');
ylabel('#count');
if ~isempty(badTrials)
    stem(h2,badTrials,sum(allBadTrialsMatrix(:,badTrials),1),'color','r');
end
subplot(h3); cla; set(h3,'nextplot','add');
stem(h3,1:length(allBadTrials),sum(allBadTrialsMatrix,2)); axis('tight'); ylabel('#count');
if ~isempty(badElecs)
    stem(h3,badElecs,sum(allBadTrialsMatrix(badElecs,:),2),'color','r');
end
view([90 -90]);

saveas(summaryFig,fullfile(folderSegment,[monkeyName expDate protocolName 'summmaryBadTrials.fig']),'fig');
saveas(summaryFig,[monkeyName expDate protocolName 'summmaryBadTrials.fig'],'fig');

%**************************************************************************
% fill the badTrials summary sheet #TODO
%**************************************************************************




end
