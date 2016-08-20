

function [X,psthVals,xs,spikeChannelNumber] = getMeanFiringRate(folderName,electrodeList,goodPos,timebin)

folderSpikes = fullfile(folderName,'segmentedData','Spikes');
folderLFP = fullfile(folderName,'segmentedData','LFP');

load(fullfile(folderSpikes,'spikeInfo.mat'));
load(fullfile(folderLFP,'lfpInfo.mat'));

if ~exist('electrodeList','var')
    electrodeList = neuralChannelsStored;
end

if ~exist('timebin','var')
    timebin = 10;
end

spikeChannelNumber = intersect(electrodeList,neuralChannelsStored); % channels in the electrodeList which have recorded spikes
n = length(spikeChannelNumber);

% initialize
load(fullfile(folderSpikes,['elec' num2str(spikeChannelNumber(1)) '_SID' num2str(SourceUnitID(1)) '.mat']));
[~,xs] = getPSTH(spikeData(goodPos),timebin,[timeVals(1) timeVals(end)]);
psthVals = zeros(n,length(xs));
X = cell(1,n);

for i=1:n
    % Get the data
    clear signal spikeData
    load(fullfile(folderSpikes,['elec' num2str(spikeChannelNumber(i)) '_SID' num2str(SourceUnitID(1)) '.mat']));
    [psthVals(i,:),xs] = getPSTH(spikeData(goodPos),timebin,[timeVals(1) timeVals(end)]);
    X{i} = spikeData(goodPos);
end

end