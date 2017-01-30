% Display the grid level waveforms
% Vinay Shirhatti, 04 Feb 2015
% *************************************************************************

function displayGridActivity(subjectName,expDate,protocolName,folderSourceString,gridType)

%========================Define folder paths===============================
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');
folderSpikeSegment = fullfile(folderSegment,'Segments');

%=========================Load Data Info===================================
% load LFP Information
[analogChannelsStored,timeVals,~,analogInputNums] = loadlfpInfo(folderLFP);
% load Spikes Information
[neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes);

[~,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums);

%========================Bad Trials========================================
% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
existsBadTrialFile = 0;
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    load(badTrialFile,'allBadTrials','badTrials','nameElec','checkTheseElectrodes');
    disp([num2str(length(badTrials)) ' bad trials']);
    existsBadTrialFile = 1;
end

ignoreTrialNum = [];
badTrials = union(badTrials,ignoreTrialNum);
disp(['Ignoring trials: ' num2str(ignoreTrialNum)]);

useAllBadTrials = 1; % set this if bad trials are to be selected as per the
% corresponding electrode and unset this if the common bad trials are to be
% used

%========================Make the Grid=====================================

% Make the grid for plotting
figure;

if strcmpi(gridType,'ECoG')
    numRows=8;numCols=10;
elseif strcmpi(gridType,'Microelectrode')
    numRows=10;numCols=10;
else
    numRows=10;numCols=11;
end

gridPos = [0.02 0.02 0.96 0.96];
gapSmall = 0.01;
hElecHandles = getPlotHandles(numRows,numCols,gridPos,gapSmall);

numElectrodes = length(analogChannelsStored);

for ne = 1:numElectrodes
    
    % Read analog channel/electrode 1
    analogChannelString = analogChannelStringArray{ne};

    [k,j] = electrodePositionOnGrid(analogChannelsStored(ne),gridType,subjectName);

    spikeChannelNumber = [];
    unitID = [];

    if strncmp(gridType,'Microelectrode',5)
        spikeChannelNumber = neuralChannelsStored(ne);
        unitID = SourceUnitIDs(ne);
    end
    
    % Get the data
    clear signal analogData
    load(fullfile(folderLFP,[analogChannelString '.mat']));
    goodPos = 1:size(analogData,1);
    
    elecIndex1 = find(checkTheseElectrodes==analogChannelsStored(ne));
    
    % Select good trials as per the electrode(s)
    if useAllBadTrials && existsBadTrialFile && ~isempty(elecIndex1)

        elecBadTrials = allBadTrials{elecIndex1};
        disp(['No. of Bad trials for ' analogChannelString ': ' num2str(length(elecBadTrials))]);

        goodPos = setdiff(goodPos,elecBadTrials);
    elseif existsBadTrialFile
        goodPos = setdiff(goodPos,badTrials);
    end
    disp([analogChannelString 'pos: ' num2str(k) ',' num2str(j) ',n=' num2str(length(goodPos))]);
    
    if isempty(elecIndex1)
        plot(hElecHandles(k,j),timeVals,analogData(goodPos,:),'r');
    else
        plot(hElecHandles(k,j),timeVals,analogData(goodPos,:),'k');
    end
        
    set(hElecHandles(k,j),'Nextplot','add');
    if useAllBadTrials && existsBadTrialFile && ~isempty(elecIndex1)
        if ~isempty(elecBadTrials)
            plot(hElecHandles(k,j),timeVals,analogData(elecBadTrials,:),'g');
        end
    elseif existsBadTrialFile
        if ~isempty(badTrials)
            plot(hElecHandles(k,j),timeVals,analogData(badTrials,:),'g');
        end
    end
    
    text(0.1,0.3,num2str(analogChannelsStored(ne)),'unit','normalized','fontsize',12,'Parent',hElecHandles(k,j));
    
end

tmin = timeVals(1); tmax = timeVals(end);
ymin = -1200; ymax = 600;

set(hElecHandles,'xlim',[tmin tmax]);
set(hElecHandles,'ylim',[ymin ymax]);

%**************************************************************************
% Plot Spikes

if strncmp(gridType,'Microelectrode',5)
    figure;

    if strcmpi(gridType,'ECoG')
        numRows=8;numCols=10;
    elseif strcmpi(gridType,'Microelectrode')
        numRows=10;numCols=10;
    else
        numRows=10;numCols=11;
    end

    gridPos = [0.02 0.02 0.96 0.96];
    gapSmall = 0.01;
    hElecHandles = getPlotHandles(numRows,numCols,gridPos,gapSmall);

    numElectrodes = length(neuralChannelsStored);

    for ne = 1:numElectrodes

        % Read analog channel/electrode 1
        analogChannelString = analogChannelStringArray{ne};

        [k,j] = electrodePositionOnGrid(analogChannelsStored(ne),gridType,subjectName);

        spikeChannelNumber = [];
        unitID = [];

        spikeChannelNumber = neuralChannelsStored(ne);
        unitID = SourceUnitIDs(ne);
        
        clear meanSpikeWaveform
        if ~isempty(spikeChannelNumber)
            load(fullfile(folderSegment,'Segments',['elec' num2str(spikeChannelNumber) '.mat']));
            plot(hElecHandles(k,j),segmentData,'k');
            meanSpikeWaveform = mean(segmentData,2);
            set(hElecHandles(k,j),'Nextplot','add');
            plot(hElecHandles(k,j),meanSpikeWaveform,'color',[0.6 0.8 0.8],'Linewidth',2);
            set(hElecHandles(k,j),'Nextplot','replace');
        end
        disp([analogChannelString 'pos: ' num2str(k) ',' num2str(j) ',n=' num2str(size(segmentData,2))]);
    end
    
    text(0.1,0.2,num2str(analogChannelsStored(ne)),'unit','normalized','fontsize',12,'Parent',hElecHandles(k,j));
    
end

set(hElecHandles,'xlim',[0 50]);

end
