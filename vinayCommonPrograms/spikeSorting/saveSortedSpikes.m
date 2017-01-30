% save sorted spikes
% Vinay Shirhatti, 17 Jan 2017
%**************************************************************************

function saveSortedSpikes(subjectName,expDate,protocolName,folderSourceString,gridType,neuralChannelsToStore,goodStimTimes,timeStartFromBaseLine,deltaT)

if isempty(neuralChannelsToStore)
    neuralChannelsToStore = 1:96;
end

overwriteSpikeInfo = 0;
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
% Get folders
folderSegment = fullfile(folderName,'segmentedData');
folderSpikes = fullfile(folderSegment,'Spikes');
folderSpikeSegment = fullfile(folderSegment,'Segments');

folderSortedSegments = fullfile(folderSpikeSegment,'sorted');
outputFolder = fullfile(folderSpikes,'sorted');
makeDirectory(outputFolder);

analysisOnsetTimes = goodStimTimes + timeStartFromBaseLine;
totalStim = length(analysisOnsetTimes);
count = 0;

for i=1:length(neuralChannelsToStore)
    clear elecsortedSegments
    elecsortedSegments = fullfile(folderSortedSegments,['elec' num2str(neuralChannelsToStore(i)) '.mat']);

    if exist(elecsortedSegments,'file')
        clear unitID unitIDSorted timeStamp
        load(elecsortedSegments);

        clear uniqueUnitIDs
        uniqueUnitIDs = unique(unitID);
        for a=1:length(uniqueUnitIDs)
            count = count+1;
            u = uniqueUnitIDs(a);
            disp([neuralChannelsToStore(i) u]);
            neuralChannelsStored(count) = neuralChannelsToStore(i);
            SourceUnitID(count) = u;

            clear spikeData
            spikeData = cell(1,totalStim);
            clear timeStampUnit
            timeStampUnit = timeStamp(unitID==u);
            for j=1:totalStim
                spikeData{j} = timeStampUnit(intersect(find(timeStampUnit>=analysisOnsetTimes(j))...
                    ,find(timeStampUnit<analysisOnsetTimes(j)+deltaT)));
                if ~isempty(spikeData{j})
                    spikeData{j} = spikeData{j} - goodStimTimes(j);
                end
            end
            save(fullfile(outputFolder,['elec' num2str(neuralChannelsToStore(i)) ... 
                '_SID' num2str(u) '.mat']),'spikeData');
        end

        if exist('unitIDSorted','var')
            clear uniqueUnitIDs
            uniqueUnitIDs = setdiff(unique(unitIDSorted),255); % ignore unitIDSorted 255, since that has already been saved with unitID
            for a=1:length(uniqueUnitIDs)      
                count = count+1;
                u = uniqueUnitIDs(a);
                disp([neuralChannelsToStore(i) u]);
                neuralChannelsStored(count) = neuralChannelsToStore(i);
                SourceUnitID(count) = u; 

                clear spikeData
                spikeData = cell(1,totalStim);
                clear timeStampUnit
                timeStampUnit = timeStamp(unitIDSorted==u);

                for j=1:totalStim
                    spikeData{j} = timeStampUnit(intersect(find(timeStampUnit>=analysisOnsetTimes(j))...
                        ,find(timeStampUnit<analysisOnsetTimes(j)+deltaT)));
                    if ~isempty(spikeData{j})
                        spikeData{j} = spikeData{j} - goodStimTimes(j);
                    end
                end
                save(fullfile(outputFolder,['elec' num2str(neuralChannelsToStore(i)) ... 
                    '_SID' num2str(u) '.mat']),'spikeData');
            end
        end

    else
        disp(['No sorted segments stored for electrode ' num2str(neuralChannelsToStore(i))]);
    end
end

spikeInfoFile = fullfile(outputFolder,'spikeInfo.mat');
if ~overwriteSpikeInfo
    if exist(spikeInfoFile,'file')
        oldFile = load(spikeInfoFile);
        neuralChannelsStored = cat(1,oldFile.neuralChannelsStored,neuralChannelsStored);
        SourceUnitID = cat(1,oldFile.SourceUnitID,SourceUnitID);
    end
end
% Write sorted Spike information
save(fullfile(outputFolder,'spikeInfo.mat'),'neuralChannelsStored','SourceUnitID');
    
end