% display the ISIs distribution and the spike train for the selected spike
% segments
% Vinay Shirhatti, 30 January 2017
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function [handles] = displayISIandSpikeTrain(folderSegments,electrode,varargin)

if sum(strcmpi('sortedspikes',varargin))
    sortedspikes=1;
else
    sortedspikes=0;
end

if sum(strcmpi('rawclusters',varargin))
    rawclusters=1;
    clusters = varargin{find(strcmp(varargin,'rawclusters'))+1};
else
    rawclusters=0;
end

if sum(strcmpi('clusterOutliers',varargin))
    clusterOutliers = varargin{find(strcmp(varargin,'clusterOutliers'))+1};
else
    clusterOutliers{1} = 0;
end

if sum(strcmpi('chooseClusters',varargin))
    chooseClusters = varargin{find(strcmp(varargin,'chooseClusters'))+1};
else
    chooseClusters = [];
end

if sum(strcmpi('handles',varargin))
    handles = varargin{find(strcmp(varargin,'handles'))+1};
    nohandles = 0;
else
    handles = [];
    nohandles = 1;
end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

if sortedspikes
    folderSpikeSegments = fullfile(folderSegments,'sorted');
else
    folderSpikeSegments = folderSegments;
end

elecsortedSegments = fullfile(folderSpikeSegments,['elec' num2str(electrode) '.mat']);

if ~(rawclusters || exist(elecsortedSegments,'file')) % this will be true only when both are 0
    disp('!!???!!! No spike segments for this unit!');
    handles = [];
    return;
    
elseif rawclusters % when clusters are passed

    elecSegments = load(fullfile(folderSegments,['elec' num2str(electrode) '.mat']));
    timeStamp = elecSegments.timeStamp;
    
    if isempty(chooseClusters)
        chooseClusters = unique(clusters);
    end
    
    if clusterOutliers{1}==0
        clear clusterOutliers
        clusterOutliers = cell(1,length(chooseClusters));
    end
    
    ignoreClusters = setdiff(unique(clusters),chooseClusters);
    if ~isempty(ignoreClusters)
        for j=1:length(ignoreClusters)
            cluster = find(clusters==ignoreClusters(j));
            unitID(cluster) = 255; % put all rejected clusters into SID255
            unitIDSorted(cluster) = 255;
            outlierID(clusterOutliers{ignoreClusters(j)}) = 255;
        end
    end
    
    for i=1:length(chooseClusters)
        cluster = find(clusters==chooseClusters(i));
        cluster = setdiff(cluster,clusterOutliers{chooseClusters(i)});
        unitID(cluster) = 0;  % mark all sorted and accepted spike clusters as SID0
        unitIDSorted(cluster) = i;
        unitID(clusterOutliers{chooseClusters(i)}) = 255; % ouliers marked as noisy spikes
        unitIDSorted(clusterOutliers{chooseClusters(i)}) = 255; % ouliers marked as noisy spikes
        outlierID(clusterOutliers{chooseClusters(i)}) = i;
    end
    
elseif exist(elecsortedSegments,'file') % read saved spike unit data
    
    clear unitID unitIDSorted timeStamp
    elecSegments = load(elecsortedSegments);
    timeStamp = elecSegments.timeStamp;
    unitID = elecSegments.unitID;
    unitIDSorted = elecSegments.unitIDSorted;
end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% plot the ISI and spike trains

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
    hfig = figure('name',[folderSegments '_elec' num2str(electrode) ' sorted clusters'],'numbertitle','off');
end    
yheightISI = (0.9/numUnits)*0.55; yheightSpkTrn = (0.9/numUnits)*0.15; ygap = (0.9/(2*numUnits-1))*0.3;  
ygapSpkTrn = yheightSpkTrn+ygap; ygapISI = yheightISI+ygap;

ucolor = gray(length(uniqueUnitIDs)+2);
ucolor = ucolor(end-2:-1:1,:);% remove the whitest colors and bring the light shades on top => black will be for the noisy, unsorted ID 255
c1=0; c2=0; count=0;
for a=1:numUnits
    count=count+1;
    u = allUnitIDs(a);
    disp([electrode u]);
    clear timeStampUnit
    if sum(u==uniqueUnitIDs) % this ID belongs to unitID
        c1 = c1+1;
        timeStampUnit = timeStamp(unitID==u);
        plotcolor = ucolor(c1,:);
    else % this belongs to unitIDSorted
        c2=c2+1;
        timeStampUnit = timeStamp(unitIDSorted==u);
        plotcolor = ucolorsorted(c2,:);
    end

    % Create handles for ISI and Spike train plots if not already
    % passed as an argument
    if nohandles
        hISI = getPlotHandles(1,1,[0.07 0.96-(a*(yheightISI)+(a-1)*(2*ygap+yheightSpkTrn)) 0.9 yheightISI]);
        hSpkTrn = getPlotHandles(1,1,[0.07 0.96-a*(yheightSpkTrn+yheightISI)-(2*a-1)*ygap 0.9 yheightSpkTrn]);
        handles.hISI(count,:) = hISI;
        handles.hSpkTrn(count,:) = hSpkTrn;
    else
        hISI = handles.hISI(count,:);
        hSpkTrn = handles.hSpkTrn(count,:);
    end

    subplot(hISI); cla; set(hISI,'nextplot','add');
    numbins = 100; binEdges = 0:0.005:0.5; % binwidth: 5ms, 100 bins from 0 to 500ms
    hhist = histogram(hISI,diff(timeStampUnit),'normalization','count','facecolor',plotcolor,'numBins',numbins,'BinEdges',binEdges);
    text(0.75,0.9,['unitID:' num2str(u)],'unit','normalized','color',plotcolor,'parent',hISI);
    text(0.9,0.9,['n=' num2str(length(timeStampUnit))],'unit','normalized','color',0.8*plotcolor,'parent',hISI);

    subplot(hSpkTrn); cla; set(hSpkTrn,'nextplot','add');
    set(hSpkTrn,'xlim',[timeStamp(1) timeStamp(end)]);
    set(hSpkTrn,'ylim',[0 1]);
    set(hSpkTrn,'yticklabelmode','manual'); set(hSpkTrn,'ytick',[]); set(hSpkTrn,'yticklabel',[]);
    for l=1:length(timeStampUnit)
        line([timeStampUnit(l) timeStampUnit(l)],[0 1],'color',plotcolor,'linewidth',1,'parent',hSpkTrn);
    end
    drawnow;
end

end