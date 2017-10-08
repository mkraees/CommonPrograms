% supporting function to read out firing rates, psth for specified
% electrodes and specifed protocol, corresponding to the specified goodPos
% Vinay Shirhatti, October 2017
% 
% modified to include option to calculate psth for required electrodes and 
% their multiple unitIDs; for eg. elec45_SID0, elec45_SID1 etc.
% 01 February 2017
% _________________________________________________________________________
% ------------------------input arguments----------------------------------
% ======= required:
% folderName: path to the folder that contains the protocol data
%
% ======= varargin:
% 'useSortedSpikesKMeans': include this to use sorted spikes (kmeans
%                           sorting in MATLAB, stored in /sorted folder
%                           separately [Vinay], not the 'spikesort' sorting
%                           , whose SIDs are stored in the original folder
%                           directly)
%
% 'removeUnitID255' : include this to use sorted spikes barring unitID255,
%                     otherwise unitID255 is used as well
%
% 'electrodeList',electrodeList : include this pair to specify the
%                                   electrodes to process
%
% 'unitIDList',unitIDList : include this pair to specify the list of 
%                           unitIDs. When both electrodeList and unitIDList
%                           are passed they must be passed as a matching
%                           list i.e. of same length and such that
%                           electrodList(n) and unitIDList(n) form the
%                           required pair
%
% 'useAllUnitIDs' : include this to process all available unitIDs for the
%                   specified electrodes. If this is not included then only
%                   unitID 0 is considered for processing
%
% 'useAsMultiunit' : include this option to combine all valid unitIDs for
%                    an electrode and treat all units together as a
%                    multiunit
%
% 'snrThreshold' : minimum acceptable SNR for the spike unit
%
% 'goodPos' : the good stimulus positions (indices) to calculate spiking 
%               data and psth for
%
% 'timebin' : width of the time bin in milliseconds for psth calculation
%
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


function [X,psthVals,xs,electrodeList,unitIDList] = getMeanFiringRate(folderName,varargin)

%__________________________________________________________________________
% Decide the path and load info
if sum(strcmpi('useSortedSpikesKMeans',varargin))
    useSortedSpikesKMeans=1;
else
    useSortedSpikesKMeans=0;
end

if ~useSortedSpikesKMeans
    folderSpikes = fullfile(folderName,'segmentedData','Spikes');
else
    folderSpikes = fullfile(folderName,'segmentedData','Spikes','sorted');
end
folderLFP = fullfile(folderName,'segmentedData','LFP');
load(fullfile(folderSpikes,'spikeInfo.mat'));
load(fullfile(folderLFP,'lfpInfo.mat'));

%__________________________________________________________________________
% read goodPos and timebin
if sum(strcmpi('goodPos',varargin))
    goodPos = varargin{find(strcmp(varargin,'goodPos'))+1};
    useAllGoodPos = 0;
else
    goodPos = 0; % just a dummy
    useAllGoodPos = 1;
end

if iscell(goodPos)
    goodPos = cell2mat(goodPos);
end

if sum(strcmpi('timebin',varargin))
    timebin = varargin{find(strcmp(varargin,'timebin'))+1};
else
    timebin = 10; % default: 10 ms
end

%__________________________________________________________________________
% check for the inputs or load defaults
% read electrodes and SourceUnitIDs
if sum(strcmpi('useAllunitIDs',varargin))
    useAllunitIDs = 1; % use all the available unitIDs
else
    useAllunitIDs = 0;
end

if sum(strcmpi('electrodeList',varargin)) % check if a list of electrodes has been passed
    electrodeList = varargin{find(strcmp(varargin,'electrodeList'))+1};
    useAllElectrodes = 0;
    
    if useAllunitIDs % use all the stored unitIDs for the specified electrode
        electrodeList = unique(electrodeList); % the user may have passed a set of electrodes with repeats, remove such repeats
        newElecList = [];
        for e=1:length(electrodeList)
            if sum(ismember(neuralChannelsStored,electrodeList(e)))
                % check the number of unitIDs stored for each electrode and
                % repeat the entries for the electrodes accordingly
                newElecList = cat(2,newElecList,repmat(electrodeList(e),1,sum(ismember(neuralChannelsStored,electrodeList(e)))));
            else
                newElecList = cat(2,newElecList,electrodeList(e));
            end
        end
        electrodeList = newElecList;
    end
else
    electrodeList = neuralChannelsStored;
    useAllElectrodes = 1;
end

if sum(strcmpi('unitIDList',varargin))
    unitIDList = varargin{find(strcmp(varargin,'unitIDList'))+1};
else
    if useAllunitIDs
        if useAllElectrodes
            unitIDList = SourceUnitID; % the entire list
        else
            unitIDList = [];
            uniqueElec = unique(electrodeList);
            for e=1:length(uniqueElec)
                if sum(ismember(neuralChannelsStored,uniqueElec(e)))
                    unitIDList = cat(2,unitIDList,SourceUnitID(ismember(neuralChannelsStored,uniqueElec(e)))); % entries corresponding to the 
                else
                    unitIDList = cat(2,unitIDList,0);
                end
            end
        end
    else % use only unit ID 0
        if ~isempty(electrodeList)
            unitIDList = zeros(1,length(electrodeList));
        else
            unitIDList = [];
        end
    end
end

% remove all unitIDs equal to 255 which contain noisy spikes, if asked to
if sum(strcmpi('removeUnitID255',varargin))
    removeUnitID255 = varargin{find(strcmp(varargin,'removeUnitID255'))+1};
else
    removeUnitID255 = 0;
end
if removeUnitID255
    electrodeList(unitIDList==255)=[];
    unitIDList(unitIDList==255)=[];
end

%****************
% check the SNR values and select units based on a minimum specified
% threshold for acceptable SNR
if sum(strcmpi('snrScreening',varargin)) % do SNR screening only if asked for
    if sum(strcmpi('snrThreshold',varargin)) % check for the snr threshold
        snrThreshold = varargin{find(strcmp(varargin,'snrThreshold'))+1};
    else
        snrThreshold = 2.5; % default minimum acceptable SNR
    end
    unitSNR = zeros(length(electrodeList),1);
    for e=1:length(electrodeList)
        clear segmentData unitID unitSegmentData
        load(fullfile(folderName,'segmentedData','Segments',['elec' num2str(electrodeList(e)) '.mat']));
        unitSegmentData = segmentData(:,unitID==unitIDList(e));
        unitSNR(e,:) = getSNR(unitSegmentData);
    end
    % remove all unitIDs with low SNR
    electrodeList(unitSNR<snrThreshold)=[];
    unitIDList(unitSNR<snrThreshold)=[];
end

%***********
% Check the mean number of spikes across units and select units based on
% minimum acceptable number
if sum(strcmpi('psthScreening',varargin)) % do psth based screening only if asked for
    maxpsthUnit = zeros(length(electrodeList),1);
    for e=1:length(electrodeList)
        clear spikeData
        load(fullfile(folderSpikes,['elec' num2str(electrodeList(e)) '_SID' num2str(unitIDList(e)) '.mat']),'spikeData');
        badTrialsFile = fullfile(folderName,'segmentedData','badTrials.mat');
        if exist(badTrialsFile,'file')
            load(badTrialsFile,'badTrials');
        else
            badTrials = [];
        end
        goodRepeats = setdiff(1:length(spikeData),badTrials);
        maxpsthUnit(e,:) = max(getPSTH(spikeData(goodRepeats),timebin,[timeVals(1) timeVals(end)]));
    end
    % remove all unitIDs with low number of spikes (psth)
    if sum(strcmpi('psthThreshold',varargin)) % check for the psth threshold
        psthThreshold = varargin{find(strcmp(varargin,'psthThreshold'))+1};
    else
        psthThreshold = 1; % default minimum acceptable psth
    end
    electrodeList(maxpsthUnit<psthThreshold)=[]; % remove the units with maxpsth<threshold
    unitIDList(maxpsthUnit<psthThreshold)=[];
end

%**** screening based on numSpikes
% if sum(strcmpi('numSpikesScreening',varargin)) % do numSpikes based screening only if asked for
%     if sum(strcmpi('spikesRange',varargin)) % check for the numSpikes threshold
%         spikesRange = varargin{find(strcmp(varargin,'spikesRange'))+1};
%     else
%         spikesRange{1} = [-0.5 0]; % check the number of spikes in this time range
%         spikesRange{2} = [0 0.25]; % check the number of spikes in this time range
%         spikesRange{3} = [0.5 0.75]; % check the number of spikes in this time range
%     end
%     unitNumSpikes = zeros(length(electrodeList),length(spikesRange));
%     for e=1:length(electrodeList)
%         clear spikeData unitNumSpikesAllPos
%         load(fullfile(folderSpikes,['elec' num2str(electrodeList(e)) '_SID' num2str(unitIDList(e)) '.mat']));
%         [~,xs] = getPSTH(spikeData,timebin,[timeVals(1) timeVals(end)]);
%         for i=1:length(spikesRange)
%             unitNumSpikes(e,i) = sum(cell2mat(cellfun(@(x) sum(x>=spikesRange{i}(1) & x<=spikesRange{i}(2)), spikeData,'uniformOutput',0)));
%         end
%     end
%     % remove all unitIDs with low number of spikes
%     if sum(strcmpi('numSpikesThreshold',varargin)) % check for the numSpikes threshold
%         numSpikesThreshold = varargin{find(strcmp(varargin,'numSpikesThreshold'))+1};
%     else
%         numSpikesThreshold = 100; % default minimum acceptable number of spikes
%     end
%     boolSpikeScreening = unitNumSpikes<numSpikesThreshold; % set the conditions where numSpikes<threshold 
%     boolSpikeScreening = prod(boolSpikeScreening,2); % set the flag for units which show numSpikes<threshold for all spikesRanges
%     electrodeList(boolSpikeScreening)=[]; % remove them
%     unitIDList(boolSpikeScreening)=[];
% end

% set the flag to convert single units to a multiunit for an electrode
if sum(strcmpi('useAsMultiunit',varargin))
    useAsMultiunit = varargin{find(strcmp(varargin,'useAsMultiunit'))+1}; % combine all valid unitIDs for an electrode as a multiunit
else
    useAsMultiunit = 0;
end

%__________________________________________________________________________
% Do the calculations now
n = length(electrodeList);
spikesExist = ismember(electrodeList,neuralChannelsStored);
% initialize
if n~=0 && sum(spikesExist) % go ahead if there is at least one electrode with stored spikes
    load(fullfile(folderSpikes,['elec' num2str(electrodeList(find(spikesExist==1,1))) '_SID' num2str(SourceUnitID(1)) '.mat']));
    [~,xs] = getPSTH(spikeData(goodPos),timebin,[timeVals(1) timeVals(end)]);
    psthVals = zeros(n,length(xs));
    if useAllGoodPos
        goodPos = size(spikeData,2);
    end
    
    if useAsMultiunit
        uniqueElec = unique(electrodeList);
        n = length(uniqueElec);
        X = cell(1,n);
        psthVals = zeros(n,length(xs));
        for i=1:n
            if sum(ismember(electrodeList(i),neuralChannelsStored))
                clear spikeDataMUA
                elecIndices = find(uniqueElec(i)==electrodeList);
                elecunitIDs = unitIDList(elecIndices);
                spikeDataMUA = cell(1,length(elecunitIDs));
                for j=1:length(elecunitIDs)
                    clear spikeData
                    disp(['elec' num2str(electrodeList(elecIndices(1))) '_SID' num2str(unitIDList(elecIndices))]);
                    load(fullfile(folderSpikes,['elec' num2str(electrodeList(elecIndices(j))) '_SID' num2str(unitIDList(elecIndices(j))) '.mat']));
                    spikeDataMUA{j} = spikeData(goodPos);
                end
                clear spikeData s1
                % concatenate spikeData for all unitIDs for an electrode
                s1=reshape([spikeDataMUA{:}],length(goodPos),length(elecunitIDs));
                spikeData = arrayfun(@(Y) vertcat(s1{Y,:}),1:length(goodPos),'un',0);
                [psthVals(i,:),xs] = getPSTH(spikeData,timebin,[timeVals(1) timeVals(end)]);
                X{i} = spikeData;
            else
                psthVals(i,:) = zeros(1,length(xs));
                X{i} = [];
            end
        end
    else % treat each unitID as a single unit
        X = cell(1,n);
        for i=1:n
            if sum(ismember(electrodeList(i),neuralChannelsStored))
                % Get the data
                clear spikeData
                disp(['elec' num2str(electrodeList(i)) '_SID' num2str(unitIDList(i))]);
                load(fullfile(folderSpikes,['elec' num2str(electrodeList(i)) '_SID' num2str(unitIDList(i)) '.mat']));
                [psthVals(i,:),xs] = getPSTH(spikeData(goodPos),timebin,[timeVals(1) timeVals(end)]);
                X{i} = spikeData(goodPos);
            else
                psthVals(i,:) = zeros(1,length(xs));
                X{i} = [];
            end
        end
    end
else
    electrodeList=0;
    load(fullfile(folderSpikes,['elec' num2str(neuralChannelsStored(1)) '_SID' num2str(SourceUnitID(1)) '.mat']));
    [~,xs] = getPSTH(spikeData(goodPos),timebin,[timeVals(1) timeVals(end)]);
    psthVals = zeros(1,length(xs));
    X{1} = 0;
end

end