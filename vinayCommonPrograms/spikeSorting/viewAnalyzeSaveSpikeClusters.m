% View, analyze and store spike clusters
% Vinay Shirhatti, 16 May 2016
%
% modified from viewAnalyzeSpikeClusters.m, Jan 2017
% This code launches a GUI which can compute spike clusters based on
% k-means clustering method, for the selected protocol and electrode
% One can inspect the clusters and sort and save new units
% *************************************************************************


function viewAnalyzeSaveSpikeClusters(folderSourceString)

poolobj = parpool('local');
%**************************************************************************
% Defaults
%**************************************************************************
numClusters = 3;
chooseClusters = 1:numClusters;
outlierCutoff = 95;
clusters = [];
centroids = [];
sumd = [];
sortInfo = [];
clusterOutliers=[];
unitID=[];
unitIDSorted=[];
timeStamp=[];
hISISpkTrn = [];
hSpikeSegments = [];
%**************************************************************************
figure;

if ~exist('folderSourceString','var')
    labUbuntu = 1;
    if isunix
        if labUbuntu
            folderSourceString = '/media/vinay/SRLHD02M/';
        else
            folderSourceString = '/media/store/';
        end
    else
        folderSourceString = 'K:\';
    end
end

gridType = 'Microelectrode';
% default initialization
subjectNames = {'alpa','kesari'};
subjectName = subjectNames{1};
[expDates,protocolNames] = eval(['allProtocols' upper(subjectName(1)) subjectName(2:end) gridType]);
expDate = expDates{1}; protocolName = protocolNames{1}; pindices = [];
protocolList = cell(1,length(protocolNames));
for i=1:length(protocolNames)
    protocolList{i} = [num2str(i) '.' expDates{i} protocolNames{i}];
end
pindices = 1:length(protocolNames);
[uniqueExpDates,ui] = unique(expDates);
uniqueExpDates = ['all' expDates(sort(ui))];
%**************************************************************************
% Display Options, panels
%**************************************************************************

% Display main options
% Gaps
xGap = 0.01; yGap = 0.01;
% Fonts
fontSizeSuperTiny = 6; fontSizeTiny = 8; fontSizeSmall = 9; 
fontSizeMedium = 10; fontSizeLarge = 12; fontSizeExtraLarge = 16;
% Background colour
backgroundColor = 'w';

% ~~~~~~~~~~~~~~~~~~~~~~~Protocol Selection panel~~~~~~~~~~~~~~~~~~~~~~~~~~
startPosX = 0.65; startPosY = 0.85; widthX = 0.3; heightY = 0.12;
hProtocolPanel = uipanel('Title','Protocol Selection','fontSize',fontSizeSmall,...
    'Unit','Normalized','Position',[startPosX startPosY widthX heightY]);
% select subjectName
textWidth = 0.2; textHeight = 0.45; xgap=0.03; ygap = 0.03;
uicontrol('Parent',hProtocolPanel,'Unit','Normalized',...
    'Position',[0.05 1-(textHeight+ygap) textWidth textHeight],...
    'Style','text','String','subject','FontSize',fontSizeSmall);
hSubjectName = uicontrol('Parent',hProtocolPanel,'Unit','Normalized',...
    'Position',[0.05+textWidth+xgap 1-(textHeight+ygap) textWidth textHeight],...
    'Style','popup','String',subjectNames,'FontSize',fontSizeSmall,'Callback',{@subjectName_Callback});
% select experiment Date
uicontrol('Parent',hProtocolPanel,'Unit','Normalized',...
    'Position',[0.05+2*(textWidth+xgap) 1-(textHeight+ygap) textWidth textHeight],...
    'Style','text','String','date','FontSize',fontSizeSmall);
hExpDate = uicontrol('Parent',hProtocolPanel,'Unit','Normalized',...
    'Position',[0.05+3*(textWidth+xgap) 1-(textHeight+ygap) textWidth textHeight],...
    'Style','popup','String',uniqueExpDates,'FontSize',fontSizeSmall,'Callback',{@expDate_Callback});
% select protocol
uicontrol('Parent',hProtocolPanel,'Unit','Normalized',...
    'Position',[0.05 1-2*(textHeight+ygap) textWidth textHeight],...
    'Style','text','String','protocol','FontSize',fontSizeSmall);
hProtocolName = uicontrol('Parent',hProtocolPanel,'Unit','Normalized',...
    'Position',[0.05+(textWidth+xgap) 1-2*(textHeight+ygap) 2*textWidth textHeight],...
    'Style','popup','String',protocolList,'FontSize',fontSizeSmall,'Callback',{@protocolName_Callback});


%~~~~~~~~~~~~~~~~~~~~~~~~Show the electrode grid~~~~~~~~~~~~~~~~~~~~~~~~~~~
gridLayoutPosition = [0.65 0.57 0.3 0.25];
numRows=10;numCols=10;
gridPlotHandle = subplot('Position',gridLayoutPosition,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[],'box','on');
axes(gridPlotHandle);
dX = 1/numCols;
dY = 1/numRows;

lineXRow = zeros(2,numRows);lineYRow = zeros(2,numRows);
for i=1:numRows
    lineXRow(:,i) = [0 1]; lineYRow(:,i) = [i*dY i*dY];
end
lineXCol = zeros(2,numCols);lineYCol = zeros(2,numCols);
for i=1:numCols
    lineXCol(:,i) = [i*dX i*dX]; lineYCol(:,i) = [0 1];
end
line(lineXRow,lineYRow,'color','k'); hold on;
line(lineXCol,lineYCol,'color','k'); 
hold off;

% load the RF data
rfDataFile = [subjectName gridType 'RFData.mat']; % cutoff = 100
electrodeList = 1:96;
if exist(rfDataFile,'file')
    highRMSElectrodes = loadRFData(rfDataFile);
    electrodeList = highRMSElectrodes;
end

[~,~,electrodeArray] = electrodePositionOnGrid(1,gridType,subjectName);
% patchColor = 'w';
elecSelHandle = cell(numRows,numCols);
selectedElec = []; % stores the selected electrode

for i=1:numRows
    textY = (numRows-i)*dY + dY/2;
    for j=1:numCols
        textX = (j-1)*dX + dX/2;
        if electrodeArray(i,j)>0
            text(textX,textY,num2str(electrodeArray(i,j)),'HorizontalAlignment','center');
        end
    end
end
set(gridPlotHandle,'XTickLabel',[],'YTickLabel',[]);

for i = 1:numRows
    for j = 1:numCols
        
        patchX = (j-1)*dX;
        patchY = (numRows - i)*dY;
        patchLocX = [patchX patchX patchX+dX patchX+dX];
        patchLocY = [patchY patchY+dY patchY+dY patchY];
        
        if sum(find(electrodeArray(i,j)==highRMSElectrodes))>0
            patchColor(i,j,:) = [0.3 0.5 0.7];
        else
            patchColor(i,j,:) = [1 1 1];
        end
        
        elecSelHandle{i,j} = patch('XData',patchLocX,'YData',patchLocY,'FaceColor',patchColor(i,j,:),...
            'FaceAlpha',0.3,'ButtonDownFcn',{@selectSingleElectrode_Callback,i,j});
        
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~Analysis options~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
startPosX = 0.65; startPosY = 0.29; widthX = 0.3; heightY = 0.25;
hOptionsPanel = uipanel('Title','Analysis Options','fontSize',fontSizeSmall,...
    'Unit','Normalized','Position',[startPosX startPosY widthX heightY]);
% clustering options
textWidth2 = 0.2; textHeight2 = 0.2; xgap=0.03; ygap = 0.03;
uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05 1-(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','text','String','number of clusters','FontSize',fontSizeSmall);
hNumClusters = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+textWidth2+xgap 1-(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','edit','String',num2str(numClusters),'FontSize',fontSizeSmall,'Callback',{@numClusters_Callback});
uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+2*(textWidth2+xgap) 1-(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','text','String','clustering method','FontSize',fontSizeSmall);
hClusterMethod = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+3*(textWidth2+xgap) 1-(textHeight2/2+ygap) textWidth2 textHeight2/2],...
    'Style','popup','String','k-means','FontSize',fontSizeSmall);
clusterDistanceMethod = {'sqeuclidean','cosine','correlation','cityblock'};
hClusterDistanceMethod = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+3*(textWidth2+xgap) 1-(textHeight2+1.5*ygap) textWidth2 textHeight2/2],...
    'Style','popup','String',clusterDistanceMethod,'FontSize',fontSizeSmall);
uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05 1-2*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','text','String','outlier cutoff (0-100)','FontSize',fontSizeSmall);
hOutlierCutoff = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+textWidth2+xgap 1-2*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','edit','String',num2str(outlierCutoff),'FontSize',fontSizeSmall,'Callback',{@outlierCutoff_Callback});
hPlotClusters = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+2*(textWidth2+xgap) 1-2*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','pushbutton','String','plot','FontSize',fontSizeSmall,'Callback',{@plotClusters_Callback});
hPlotOption = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+3*(textWidth2+xgap) 1-2*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','popup','String','show all|with outliers|w/o outliers','FontSize',fontSizeSmall);
uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05 1-3*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','text','String','choose clusters','FontSize',fontSizeSmall);
hChooseClusters = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+textWidth2+xgap 1-3*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','edit','String',num2str(chooseClusters),'FontSize',fontSizeSmall,'Callback',{@chooseClusters_Callback});
hSaveClusters = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+2*(textWidth2+xgap) 1-3*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','pushbutton','String','save','FontSize',fontSizeSmall,'Callback',{@saveClusters_Callback});
hDeleteClusters = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+3*(textWidth2+xgap) 1-3*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','pushbutton','String','delete','FontSize',fontSizeSmall,'Callback',{@deleteClusters_Callback});
hClearPlot = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05 1-4*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','pushbutton','String','clear plot','FontSize',fontSizeSmall,'Callback',{@cla_Callback});
hPlotSavedSpikeSegments = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+textWidth2+xgap 1-4*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','pushbutton','String','plot saved','FontSize',fontSizeSmall,'Callback',{@plotSavedSegments_Callback});
hSegmentsType = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.05+2*(textWidth2+xgap) 1-4*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','popup','String','sorted|unsorted','FontSize',fontSizeSmall);
hExit = uicontrol('Parent',hOptionsPanel,'Unit','Normalized',...
    'Position',[0.95-textWidth2 1-4*(textHeight2+ygap) textWidth2 textHeight2],...
    'Style','pushbutton','String','exit pool','FontSize',fontSizeSmall,'Callback',{@exitpool_Callback});


%~~~~~~~~~~~~~~~~~~~~~~~~~~~Spikes Cluster plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clustersPlotPos = [0.05 0.05 0.55 0.9];
gapX = 0.04; gapY = 0.04;
hClustersPlot = getPlotHandles(ceil((numClusters+1)/2),2,clustersPlotPos,gapX,gapY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CALLBACK FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selectSingleElectrode_Callback(~,~,row,col)
    
    [~,~,electrodeArray] = electrodePositionOnGrid(1,gridType,subjectName);
    
    for i = 1:numRows
        for j = 1:numCols
            if sum(find(electrodeArray(i,j)==highRMSElectrodes))>0
                patchColor(i,j,:) = [0.3 0.5 0.7];
            else
                patchColor(i,j,:) = [1 1 1];
            end
        end
    end
    
    thisSelElec = electrodeArray(row,col);
    if isempty(selectedElec)
        selectedElec = [];
    end
    if thisSelElec == 0
        disp('Select valid Electrode');
    else
        if sum(find(selectedElec==thisSelElec))==0 % electrode not selected yet
            selectedElec = thisSelElec;
            disp(['Electrode selected: ' num2str(thisSelElec)]);
            for r=1:numRows
                for c=1:numCols
                    if r==row && c==col
                        set(elecSelHandle{row,col},'FaceColor','r');
                    else
                        set(elecSelHandle{r,c},'FaceColor',patchColor(r,c,:));
                    end
                end
            end
        else
            disp(['Deselecting Electrode: ' num2str(thisSelElec)]);
            selectedElec(selectedElec==thisSelElec) = [];
            set(elecSelHandle{row,col},'FaceColor',patchColor(row,col,:));
        end
    end
end

%**************************************************************************

function subjectName_Callback(~,~)
    
    subjectName = subjectNames{get(hSubjectName,'val')};

    if strcmpi(subjectName,'alpa')
        [expDates,protocolNames] = allProtocolsAlpaMicroelectrode;
    elseif strcmpi(subjectName,'kesari')
    %     [expDates,protocolNames,stimType] = allProtocolsKesariMicroelectrode;
        [expDates,protocolNames] = getAllProtocols(subjectName,gridType);
    else
        [expDates,protocolNames] = eval(['allProtocols' upper(subjectName(1)) subjectName(2:end) gridType]);
    end
    protocolList = cell(1,length(protocolNames));
    for i=1:length(protocolNames)
        protocolList{i} = [num2str(i) '.' expDates{i} protocolNames{i}];
    end
    pindices = 1:length(protocolNames);
    [uniqueExpDates,ui] = unique(expDates);
    uniqueExpDates = ['all' expDates(sort(ui))];
    
    hExpDate = uicontrol('Parent',hProtocolPanel,'Unit','Normalized',...
        'Position',[0.05+3*(textWidth+xgap) 1-(textHeight+ygap) textWidth textHeight],...
        'Style','popup','String',uniqueExpDates,'FontSize',fontSizeSmall,'Callback',{@expDate_Callback});
    
    hProtocolName = uicontrol('Parent',hProtocolPanel,'Unit','Normalized',...
        'Position',[0.05+(textWidth+xgap) 1-2*(textHeight+ygap) 2*textWidth textHeight],...
        'Style','popup','String',protocolList,'FontSize',fontSizeSmall,'Callback',{@protocolName_Callback});
    
    clear elecSelHandle
    elecSelHandle = cell(numRows,numCols);
    numRows = 10; numCols = 10;

    axes(gridPlotHandle);
    cla;

    dX = 1/numCols;
    dY = 1/numRows;

    lineXRow = zeros(2,numRows);lineYRow = zeros(2,numRows);
    for i=1:numRows
        lineXRow(:,i) = [0 1]; lineYRow(:,i) = [i*dY i*dY];
    end
    lineXCol = zeros(2,numCols);lineYCol = zeros(2,numCols);
    for i=1:numCols
        lineXCol(:,i) = [i*dX i*dX]; lineYCol(:,i) = [0 1];
    end
    line(lineXRow,lineYRow,'color','k'); hold on;
    line(lineXCol,lineYCol,'color','k'); 
    hold off;

    % load the RF data
    rfDataFile = [subjectName gridType 'RFData.mat']; % cutoff = 100
    electrodeList = 1:96;
    if exist(rfDataFile,'file')
        highRMSElectrodes = loadRFData(rfDataFile);
        electrodeList = highRMSElectrodes;
    end

    [~,~,electrodeArray] = electrodePositionOnGrid(1,gridType,subjectName);
    % patchColor = 'w';

    for i=1:numRows
        textY = (numRows-i)*dY + dY/2;
        for j=1:numCols
            textX = (j-1)*dX + dX/2;
            if electrodeArray(i,j)>0
                text(textX,textY,num2str(electrodeArray(i,j)),'HorizontalAlignment','center');
            end
        end
    end
    set(gridPlotHandle,'XTickLabel',[],'YTickLabel',[]);

    for i = 1:numRows
        for j = 1:numCols
            
            patchX = (j-1)*dX;
            patchY = (numRows - i)*dY;
            patchLocX = [patchX patchX patchX+dX patchX+dX];
            patchLocY = [patchY patchY+dY patchY+dY patchY];

            if sum(find(electrodeArray(i,j)==highRMSElectrodes))>0
                patchColor(i,j,:) = [0.3 0.5 0.7];
            else
                patchColor(i,j,:) = [1 1 1];
            end

            elecSelHandle{i,j} = patch('XData',patchLocX,'YData',patchLocY,'FaceColor',patchColor(i,j,:),...
                'FaceAlpha',0.3,'ButtonDownFcn',{@selectSingleElectrode_Callback,i,j});

        end
    end
    
end

%**************************************************************************

function expDate_Callback(~,~)    
    expDate = uniqueExpDates{get(hExpDate,'val')};
    pindices = 1:length(protocolNames);
    if get(hExpDate,'val')~=1
        pindices = pindices(ismember(expDates,expDate));
    else
        pindices = 1:length(protocolNames);
    end
    
    protocols = {protocolList{pindices}};
    hProtocolName = uicontrol('Parent',hProtocolPanel,'Unit','Normalized',...
        'Position',[0.05+(textWidth+xgap) 1-2*(textHeight+ygap) 2*textWidth textHeight],...
        'Style','popup','String',protocols,'FontSize',fontSizeSmall,'Callback',{@protocolName_Callback});
end

%**************************************************************************

function protocolName_Callback(~,~)
    
    pindices = pindices(get(hProtocolName,'val'));
    protocolName = protocolNames{pindices};
    
end

%**************************************************************************

function numClusters_Callback(~,~)   
    numClusters = str2double(get(hNumClusters,'String'));
    hClustersPlot = getPlotHandles(ceil((numClusters+1)/2),2,clustersPlotPos,gapX,gapY);
end

%**************************************************************************

function outlierCutoff_Callback(~,~)
    outlierCutoff = str2double(get(hOutlierCutoff,'String'));
end

%**************************************************************************

function chooseClusters_Callback(~,~)
    chooseClusters = str2double(get(hChooseClusters,'String'));
    
    folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
    % Get folders
    folderSegment = fullfile(folderName,'segmentedData');
    folderSpikeSegment = fullfile(folderSegment,'Segments');
    hISISpkTrn = displayISIandSpikeTrain(folderSpikeSegment,selectedElec,'rawclusters',clusters,'clusterOutliers',clusterOutliers,'chooseClusters',chooseClusters);
end

%**************************************************************************

function plotClusters_Callback(~,~)
    numClusters = str2double(get(hNumClusters,'String'));
    outlierCutoff = str2double(get(hOutlierCutoff,'String'));
    plotOption = get(hPlotOption,'val');
    distancemethod = clusterDistanceMethod{get(hClusterDistanceMethod,'val')};
    
    [clusters,centroids,sumd,sortInfo,clusterOutliers] = plotSpikeClusters(hClustersPlot,subjectName,expDate,protocolName,folderSourceString,gridType,selectedElec,numClusters,outlierCutoff,plotOption,distancemethod);
    
    %======================================================================
    % plot ISI and spike train
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
    % Get folders
    folderSegment = fullfile(folderName,'segmentedData');
    folderSpikeSegment = fullfile(folderSegment,'Segments');
    
    hISISpkTrn = displayISIandSpikeTrain(folderSpikeSegment,selectedElec,'rawclusters',clusters,'clusterOutliers',clusterOutliers);
% 
%     elecSegments = load(fullfile(folderSpikeSegment,['elec' num2str(selectedElec) '.mat']));
%     timeStamp = elecSegments.timeStamp;
%     
%     chooseClusters = unique(clusters);
%     
%     for i=1:length(chooseClusters)
%         cluster = find(clusters==chooseClusters(i));
%         cluster = setdiff(cluster,clusterOutliers{chooseClusters(i)});
%         unitID(cluster) = 0;  % mark all sorted and accepted spike clusters as SID0
%         unitIDSorted(cluster) = i;
%         unitID(clusterOutliers{chooseClusters(i)}) = 255; % ouliers marked as noisy spikes
%         unitIDSorted(clusterOutliers{chooseClusters(i)}) = 255; % ouliers marked as noisy spikes
%         outlierID(clusterOutliers{chooseClusters(i)}) = i;
%     end
%     
%     uniqueUnitIDs = unique(unitID);
%     if exist('unitIDSorted','var')
%         uniqueUnitIDsSorted = unique(unitIDSorted);
%         ucolorsorted = hsv(length(uniqueUnitIDsSorted)-1); % no color for 255 since it is accounted for by unitID
%     else
%         uniqueUnitIDsSorted = [];
%     end
%     
%     allUnitIDs = unique([uniqueUnitIDs uniqueUnitIDsSorted]);
%     numUnits = length(allUnitIDs);
%     
%     hfig = figure('name',[subjectName expDate protocolName '_elec' num2str(selectedElec) ' sorted clusters'],'numbertitle','off');
%     yheightISI = (0.9/numUnits)*0.55; yheightSpkTrn = (0.9/numUnits)*0.15; ygap = (0.9/(2*numUnits-1))*0.3;
%     
%     ucolor = gray(length(uniqueUnitIDs)+2);
%     ucolor = ucolor(end-2:-1:1,:);% remove the whitest colors and bring the light shades on top => black will be for the noisy, unsorted ID 255
%     c1=0; c2=0; count=0;
%     for a=1:numUnits
%         count=count+1;
%         u = allUnitIDs(a);
%         disp([selectedElec u]);
%         clear timeStampUnit
%         if sum(u==uniqueUnitIDs) % this ID belongs to unitID
%             c1 = c1+1;
%             timeStampUnit = timeStamp(unitID==u);
%             plotcolor = ucolor(c1,:);
%         else
%             c2=c2+1;
%             timeStampUnit = timeStamp(unitIDSorted==u);
%             plotcolor = ucolorsorted(c2,:);
%         end
%         
%         hISI = getPlotHandles(1,1,[0.07 0.96-(a*(yheightISI)+(a-1)*(2*ygap+yheightSpkTrn)) 0.9 yheightISI]);
%         hSpkTrn = getPlotHandles(1,1,[0.07 0.96-a*(yheightSpkTrn+yheightISI)-(2*a-1)*ygap 0.9 yheightSpkTrn]);
%         
%         subplot(hISI); cla; set(hISI,'nextplot','add');
%         hhist = histogram(hISI,diff(timeStampUnit),'normalization','count','facecolor',plotcolor,'numBins',100);
%         
%         subplot(hSpkTrn); cla; set(hSpkTrn,'nextplot','add');
%         set(hSpkTrn,'xlim',[timeStamp(1) timeStamp(end)]);
%         set(hSpkTrn,'ylim',[0 1]);
%         set(hSpkTrn,'yticklabelmode','manual'); set(hSpkTrn,'ytick',[]); set(hSpkTrn,'yticklabel',[]);
%         for l=1:length(timeStampUnit)
%             line([timeStampUnit(l) timeStampUnit(l)],[0 1],'color',plotcolor,'linewidth',1,'parent',hSpkTrn);
%         end
% 
%     end
    
end

%**************************************************************************

function saveClusters_Callback(~,~)
    chooseClusters = str2double(get(hChooseClusters,'String'));
    ignoreClusters = setdiff(1:numClusters,chooseClusters);
    
    folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
    % Get folders
    folderSegment = fullfile(folderName,'segmentedData');
    folderSpikeSegment = fullfile(folderSegment,'Segments');
    
    pathname = fullfile(folderSpikeSegment,'sorted');
    makeDirectory(pathname);
    
    elecSegments = load(fullfile(folderSegment,'Segments',['elec' num2str(selectedElec) '.mat']));
    
    outlierID = 256*ones(1,length(clusters)); % initialize the outlier IDs as 256 i.e. not belonging to any cluster
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
    savefile = fullfile(pathname,['elec' num2str(selectedElec) '.mat']);
    if exist(savefile,'file')
        choice = questdlg('Sorted segments file already exists for this electrode. Do you want to replace it?', ...
            'File alert', ...
            'Yes','No','I am not sure','No');
        % Handle response
        switch choice
            case 'Yes'
                disp(['>.>.>.>.>...Saving sorted segments...  elec' num2str(selectedElec)]);
                segmentInfo = elecSegments.segmentInfo;
                segmentData = elecSegments.segmentData;
                timeStamp = elecSegments.timeStamp;
                sampleCount = elecSegments.sampleCount;
                sortInfo.clusters = clusters;
                sortInfo.numClusters = numClusters;
                sortInfo.chooseClusters = chooseClusters;
                sortInfo.centroids = centroids;
                sortInfo.sumd = sumd;
                sortInfo.clusterOutliers = clusterOutliers;
                sortInfo.outlierCutoff = outlierCutoff;
                save(savefile,'segmentInfo','timeStamp','segmentData','sampleCount','unitID','unitIDSorted','sortInfo','outlierID');
            case 'No'
                msgbox('Ok, take your time...Not saving the new file for now');
            case 'I am not sure'
                msgbox('Ok, take your time. Go take a break if you are tired, like a stroll or a coffee or some music or a nap ...Not saving the new file for now');
        end
    else
        disp(['>.>.>.>.>...Saving sorted segments...  elec' num2str(selectedElec)]);
        segmentInfo = elecSegments.segmentInfo;
        segmentData = elecSegments.segmentData;
        timeStamp = elecSegments.timeStamp;
        sampleCount = elecSegments.sampleCount;
        sortInfo.clusters = clusters;
        sortInfo.numClusters = numClusters;
        sortInfo.chooseClusters = chooseClusters;
        sortInfo.centroids = centroids;
        sortInfo.sumd = sumd;
        sortInfo.clusterOutliers = clusterOutliers;
        sortInfo.outlierCutoff = outlierCutoff;
        save(savefile,'segmentInfo','timeStamp','segmentData','sampleCount','unitID','unitIDSorted','sortInfo','outlierID');        
    end
    
    infofileName = fullfile(pathname,'segmentInfo.mat');
    if exist(infofileName,'file')
        load(infofileName);
        
        % check if sorted segments for the selected electrode are already
        % stored and clear them if they are. New sorted segments will be
        % stored now
        if ~isempty(intersect(selectedElec,segmentChannelsStored))
            elecIndex = find(segmentChannelsStored==selectedElec);
            segmentChannelsStored(elecIndex) = [];
            numItems(elecIndex) = [];
            
            if isempty(segmentChannelsStored)
                segmentChannelsStored=[];
                numItems=[];
            end
        end
        
    else
        segmentChannelsStored=[];
        numItems=[];
    end
    uniqueUnitIDs = unique(unitID);
    if exist('unitIDSorted','var')
        uniqueUnitIDsSorted = unique(unitIDSorted);
    else
        uniqueUnitIDsSorted = [];
    end

    allUnitIDs = unique([uniqueUnitIDs uniqueUnitIDsSorted]);
    numUnits = length(allUnitIDs);

    newSegmentChannels = repmat(selectedElec,1,numUnits);
    segmentChannelsStored = cat(1,segmentChannelsStored,newSegmentChannels);

    for a=1:numUnits
        u = allUnitIDs(a);
        if sum(u==uniqueUnitIDs) % this ID belongs to unitID
            newnumItems = numel(find(unitID==u));
            numItems = cat(1,numItems,newnumItems);
        else % this belongs to unitIDSorted
            newnumItems = numel(find(unitIDSorted==u));
            numItems = cat(1,numItems,newnumItems);
        end         
    end

    save(infofileName,'segmentChannelsStored','numItems');

end

%**************************************************************************

function deleteClusters_Callback(~,~)
    chooseClusters = str2double(get(hChooseClusters,'String'));
    
    pathname = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','Spikes','sorted');
    filename = ['elec' num2str(selectedElec) '_SID'];
    
    for i=1:length(chooseClusters)
        savefile = fullfile(pathname,[filename num2str(i) '.mat']);
        if exist(savefile,'file')
            delete(savefile);
        end
    end

end

%**************************************************************************

function plotSavedSegments_Callback(~,~)
    segmentsType = get(hSegmentsType,'val');
    if segmentsType==1
        hSpikeSegments = displaySpikeSegments(subjectName,expDate,protocolName,folderSourceString,gridType,selectedElec,'sortedspikes');
    elseif segmentsType==2
        hSpikeSegments = displaySpikeSegments(subjectName,expDate,protocolName,folderSourceString,gridType,selectedElec);
    end
    
end

%**************************************************************************

function cla_Callback(~,~)
        
    claGivenPlotHandle(hClustersPlot);

    function claGivenPlotHandle(plotHandles)
        [numRows,numCols] = size(plotHandles);
        for i=1:numRows
            for j=1:numCols
                cla(plotHandles(i,j));
            end
        end
    end
end

%**************************************************************************

function exitpool_Callback(~,~)
    delete(poolobj);
end

%**************************************************************************

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SUPPORTING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function highRMSElectrodes = loadRFData(filename)
load(fullfile(filename));
end

%**************************************************************************

function [clusters,centroids,sumd,sortInfo,clusterOutliers] = plotSpikeClusters(hClustersPlot,subjectName,expDate,protocolName,folderSourceString,gridType,electrodeNum,numClusters,outlierCutoff,plotOption,distancemethod)

clear clusters centroids sumd sortInfo clusterOutliers

if ~exist('outlierCutoff','var')
    outlierCutoff = 100; % no outliers calculation
end

if ~exist('plotOption','var')
    plotOption = 1;
end

%========================Define folder paths===============================
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
% Get folders
folderSegment = fullfile(folderName,'segmentedData');
folderSpikes = fullfile(folderSegment,'Spikes');
folderSpikeSegment = fullfile(folderSegment,'Segments');

%=========================Load Data Info===================================
% load Spikes Information
% [neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes);
[segmentChannelsStored,numItems] = loadsegmentInfo(folderSpikeSegment);

if isempty(intersect(electrodeNum,segmentChannelsStored))
    disp(['Electrode ' num2str(electrodeNum) ' does not have recorded spikes']);
    clusters = [];centroids=[];sumd=[];sortInfo=[];clusterOutliers=[];
    return;
end

spikeChannelNumber = electrodeNum;
colorList = getColors(numClusters);
if plotOption==1
    showOutliers = 0; showWithoutOutliers = 0;
elseif plotOption==2
    showOutliers = 1; showWithoutOutliers = 0;
else
    showOutliers = 0; showWithoutOutliers = 1;
end

clear meanSpikeWaveform clusters
set(hClustersPlot,'Nextplot','replace');
if ~isempty(spikeChannelNumber)
    load(fullfile(folderSegment,'Segments',['elec' num2str(spikeChannelNumber) '.mat']));
    opts = statset('Display','final','UseParallel',1,'MaxIter',1000);
    if ~exist('distancemethod','var')
        distancemethod = 'sqeuclidean';
%         distancemethod = 'cosine';
    end
    numreplicates = 3;
    [clusters,centroids,sumd,distanceToCentroid] = kmeans(segmentData',numClusters,'Distance',distancemethod,'Replicates',numreplicates,'Options',opts);
    
    clusterOutliers = cell(1,numClusters);
    clusterOutliersRemoved = cell(1,numClusters);
    for ii=1:numClusters
        
        cind = find(clusters==ii);
        cdist = distanceToCentroid(cind,ii);
        clusterOutliers{ii} = cind(cdist>prctile(cdist,outlierCutoff)); % flag points with distances above a particular percentile as outliers
        % for eg: if cind(cdist>prctile(cdist,95)) then distances with the
        % highest 5% of values are outliers
        clusterOutliersRemoved{ii} = setdiff(cind,clusterOutliers{ii});
        
        prow = floor(ii/2)+1; pcol = mod(ii,2)+1;
        
        set(hClustersPlot(prow,pcol),'Nextplot','replace');
        if showOutliers
            plot(hClustersPlot(prow,pcol),segmentData(:,clusterOutliers{ii}),'color',0.2*colorList(ii,:),'Linewidth',1);
            set(hClustersPlot(prow,pcol),'Nextplot','add');
            plot(hClustersPlot(prow,pcol),segmentData(:,clusterOutliersRemoved{ii}),'color',colorList(ii,:),'Linewidth',1);
            plot(hClustersPlot(prow,pcol),mean(segmentData(:,clusterOutliersRemoved{ii}),2),'color',0.5.*colorList(ii,:),'Linewidth',2.5);
        elseif ~showOutliers && showWithoutOutliers
            plot(hClustersPlot(prow,pcol),segmentData(:,clusterOutliersRemoved{ii}),'color',colorList(ii,:),'Linewidth',1);
            set(hClustersPlot(prow,pcol),'Nextplot','add');
            plot(hClustersPlot(prow,pcol),mean(segmentData(:,clusterOutliersRemoved{ii}),2),'color',0.5.*colorList(ii,:),'Linewidth',2.5);
        else
            plot(hClustersPlot(prow,pcol),segmentData(:,clusters==ii),'color',colorList(ii,:),'Linewidth',1);
            set(hClustersPlot(prow,pcol),'Nextplot','add');
            plot(hClustersPlot(prow,pcol),mean(segmentData(:,clusters==ii),2),'color',0.5.*colorList(ii,:),'Linewidth',2.5);
        end
        text(0.1, 0.95, ['SNR: ' num2str(getSNR(segmentData(:,clusters==ii)))],'fontweight','bold','color', [0.3 0.7 0.1],'unit','normalized','parent',hClustersPlot(prow,pcol));
        text(0.8, 0.95, ['clust: ' num2str(ii)],'fontweight','bold','color', [0.9 0.2 0.3],'unit','normalized','parent',hClustersPlot(prow,pcol));
        text(0.4, 0.95, ['w/ooSNR: ' num2str(getSNR(segmentData(:,clusterOutliersRemoved{ii})))],'fontweight','bold','color', 0.8.*[0.3 0.7 0.1],'unit','normalized','parent',hClustersPlot(prow,pcol));
        sortInfo.snr(ii) = getSNR(segmentData(:,clusters==ii));
        sortInfo.snrOutliersRemoved(ii) = getSNR(segmentData(:,clusterOutliersRemoved{ii}));
        
        if strcmpi(distancemethod,'sqeuclidean')
            plot(hClustersPlot(1,1),centroids(ii,:),'color',colorList(ii,:),'Linewidth',2);
        elseif showOutliers || showWithoutOutliers
            plot(hClustersPlot(1,1),mean(segmentData(:,clusterOutliersRemoved{ii}),2),'color',0.5.*colorList(ii,:),'Linewidth',2.5);
        else
            plot(hClustersPlot(1,1),mean(segmentData(:,clusters==ii),2),'color',0.5.*colorList(ii,:),'Linewidth',2.5);
        end
            
        set(hClustersPlot(1,1),'Nextplot','add'); set(hClustersPlot(1,1),'ylim',[-70 70]);
        
        drawnow;
    end
end
disp([num2str(spikeChannelNumber) ',n=' num2str(size(segmentData,2))]);
sortInfo.distancemethod = distancemethod;
sortInfo.numreplicates = numreplicates;
sortInfo.outlierCutoff = outlierCutoff;

end

%**************************************************************************

function colorList = getColors(numColors)
%     colorList = jet(numColors);
%     colorList = copper(numColors);
    colorList = hsv(numColors);
end
