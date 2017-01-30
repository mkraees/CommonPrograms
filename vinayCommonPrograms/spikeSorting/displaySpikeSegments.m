% Plot spike clusters
% Modified from an internal function in viewAnalyzeSpikeClustersV2.m:
% plotSpikeClusters.m
% Vinay Shirhatti, 30 January 2017
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function [hSpikeSegments] = displaySpikeSegments(subjectName,expDate,protocolName,folderSourceString,gridType,electrode,varargin)

if sum(strcmpi('sortedspikes',varargin))
    sortedspikes=1;
else
    sortedspikes=0;
end

if sum(strcmpi('handles',varargin))
    hSpikeSegments = varargin{find(strcmp(varargin,'handles'))+1};
    nohandles = 0;
else
    hSpikeSegments = [];
    nohandles = 1;
end

%========================Define folder paths===============================
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
% Get folders
folderSegment = fullfile(folderName,'segmentedData');
if sortedspikes
    folderSpikeSegment = fullfile(folderSegment,'Segments','sorted');
else
    folderSpikeSegment = fullfile(folderSegment,'Segments');
end

%=========================Load Data Info===================================
% load Spikes Information
% [neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes);
[segmentChannelsStored,numItems] = loadsegmentInfo(folderSpikeSegment);

if isempty(intersect(electrode,segmentChannelsStored))
    disp(['Electrode ' num2str(electrode) ' does not have recorded spikes']);
    return;
end

%==========================================================================

clear meanSpikeWaveform clusters
if ~isempty(electrode)
    load(fullfile(folderSpikeSegment,['elec' num2str(electrode) '.mat']));
    
    uniqueUnitIDs = unique(unitID);
    if exist('unitIDSorted','var')
        uniqueUnitIDsSorted = unique(unitIDSorted);
        ucolorsorted = hsv(length(uniqueUnitIDsSorted)-1); % no color for 255 since it is accounted for by unitID
    else
        uniqueUnitIDsSorted = [];
    end

    allUnitIDs = unique([uniqueUnitIDs uniqueUnitIDsSorted]);
    numUnits = length(allUnitIDs);

    if nohandles
        hfig = figure('name',[folderSpikeSegment '_elec' num2str(electrode) ' sorted segments'],'numbertitle','off');
        % create handles
        gapX = 0.04; gapY = 0.04;
        hSpikeSegments = getPlotHandles(ceil((numUnits+1)/2),2,[0.08 0.08 0.86 0.86],gapX,gapY);
    end    
    
    ucolor = gray(length(uniqueUnitIDs)+2);
    ucolor = ucolor(end-2:-1:1,:);% remove the whitest colors and bring the light shades on top => black will be for the noisy, unsorted ID 255
    c1=0; c2=0; count=0;
    for a=1:numUnits
        count=count+1;
        u = allUnitIDs(a);
        disp([electrode u]);
        clear segmentDataUnit
        if sum(u==uniqueUnitIDs) % this ID belongs to unitID
            c1 = c1+1;
            segmentDataUnit = segmentData(:,unitID==u);
            plotcolor = ucolor(c1,:);
        else % this belongs to unitIDSorted
            c2=c2+1;
            segmentDataUnit = segmentData(:,unitIDSorted==u);
            plotcolor = ucolorsorted(c2,:);
        end
        
        prow = floor(a/2)+1; pcol = mod(a,2)+1;
        subplot(hSpikeSegments(prow,pcol)); cla; set(hSpikeSegments(prow,pcol),'nextplot','add');
        plot(hSpikeSegments(prow,pcol),segmentDataUnit,'color',plotcolor,'Linewidth',1);
        plot(hSpikeSegments(prow,pcol),mean(segmentDataUnit,2),'color',0.5.*plotcolor,'Linewidth',2.5);
        set(hSpikeSegments(prow,pcol),'yticklabelmode','auto');
        
        text(0.1, 0.95, ['SNR: ' num2str(getSNR(segmentDataUnit))],'fontweight','bold','color', [0.3 0.7 0.1],'unit','normalized','parent',hSpikeSegments(prow,pcol));
        text(0.8, 0.95, ['clust: ' num2str(u)],'fontweight','bold','color', [0.9 0.2 0.3],'unit','normalized','parent',hSpikeSegments(prow,pcol));
        text(0.8, 0.85, ['n=' num2str(size(segmentDataUnit,2))],'fontweight','bold','color', [0.5 0.2 0.9],'unit','normalized','parent',hSpikeSegments(prow,pcol));
        
        plot(hSpikeSegments(1,1),mean(segmentDataUnit,2),'color',plotcolor,'Linewidth',2.5);          
        set(hSpikeSegments(1,1),'Nextplot','add'); set(hSpikeSegments(1,1),'ylim',[-70 70]);
        
        drawnow;
    end   
end
disp([num2str(electrode) ',n=' num2str(size(segmentData,2))]);

end
