% displayOrientationMap plots the orientation map of the recording grid.

% 15 Oct 2015
% This program expects [subjectName gridType 'RFData.mat'], which is
% obtained by running electrodeSelection.

% We assume that these protocols have fixed azimuth, elevation and
% temporalFrequency. The spatialFrequency and size can take more than one
% value, but in that case the positions (fPos and sPos) must be
% specified. Default positions = 1.

% Orientation Preferences are used only for electrodes/sessions for which
% the max firing rates are higher than frCutoff and orientation selectivity
% is higher than osCutoff.

function [finalElectrodeList,finalOrientationPref,finalOrientationSelectivity,finalFiringRates] = displayOrientationMap(subjectName,expDates,protocolNames,folderSourceString,gridType,gridLims,aPos,ePos,sPos,fPos,frCutoff,osCutoff,stdPOCutoff)

if ~exist('aPos','var');                aPos=1;                         end
if ~exist('ePos','var');                ePos=1;                         end
if ~exist('sPos','var');                sPos=1;                         end
if ~exist('fPos','var');                fPos=1;                         end
if ~exist('frCutoff','var');            frCutoff = 10;                  end
if ~exist('osCutoff','var');            osCutoff = 0.1;                 end
if ~exist('stdPOCutoff','var');         stdPOCutoff = 0.5;              end 

hOSHistAllElectrodes     = subplot('Position',[0.025 0.65 0.125 0.25]);
hPOHistAllElectrodes     = subplot('Position',[0.025 0.35 0.125 0.25]);
hStdPOHistAllElectrodes  = subplot('Position',[0.025 0.05 0.125 0.25]);

hNumElectrodes           = subplot('Position',[0.2 0.65 0.3 0.25]);
hMeanOriAllElectrodes    = subplot('Position',[0.2 0.35 0.3 0.25]);
hStdPOAllElectrodes      = subplot('Position',[0.2 0.05 0.3 0.25]);

hGridPlot                = subplot('Position',[0.525 0.55 0.2 0.4],'XTickLabel',[],'YTickLabel',[]);
hRFPlotAllElectrodes     = subplot('Position',[0.75 0.55 0.2 0.4]);

%hColorMapPlot            = subplot('Position',[0.825 0.65 0.125 0.25]);
hPORFPlotAllElectrodes  = subplot('Position',[0.525 0.05 0.2 0.4]);
hOSRFPlotAllElectrodes  = subplot('Position',[0.75 0.05 0.2 0.4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show RF centers
rfData = load([subjectName gridType 'RFData.mat']);
electrodeList = rfData.highRMSElectrodes;

displaySavedRFcenters(hGridPlot,hRFPlotAllElectrodes,subjectName,gridType,electrodeList,gridLims);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numDays = length(expDates);
numElectrodes = length(electrodeList);

for i=1:numDays
    clear expDate protocolName
    [aValsUniqueAll{i},eValsUniqueAll{i},sValsUniqueAll{i},fValsUniqueAll{i},oValsUniqueAll{i}] = loadParameterCombinations(subjectName,expDates{i},protocolNames{i},folderSourceString,gridType);
    aziCenter(i)=aValsUniqueAll{i}(aPos);
    eleCenter(i)=eValsUniqueAll{i}(ePos);
end

% Check whether the values are the same across days
sizeDeg = sValsUniqueAll{1}(sPos);
spatialFreqCPD = fValsUniqueAll{1}(fPos);
oValsUnique = oValsUniqueAll{1};

for i=2:numDays
    if (sizeDeg ~= sValsUniqueAll{i}(sPos))
        error('Sizes do not match.....');
    end
    if (spatialFreqCPD ~= fValsUniqueAll{i}(fPos))
        error('Spatial Frequencies do not match.....');
    end
    if (oValsUnique ~= oValsUniqueAll{i})
        error('Orientations do not match.....');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orientation selectivity (OS) is computed in two cases. If the stimulus
% size is small, we need to find, for each electrode, only the days for
% which the stimulus was near the center of the RF and estimate the OS for
% those days only. If, however, the stimulus size is large, we compute the OS
% for all days.

sizeCutoffDeg = 2; % size cutoff
dRange = 0.2; % distance range

for i=1:numElectrodes
    clear d a e
    a = rfData.rfStats(electrodeList(i)).meanAzi;
    e = rfData.rfStats(electrodeList(i)).meanEle;
    
    for j=1:numDays
        d(j) = sqrt((a-aziCenter(j))^2 + (e-eleCenter(j))^2);
    end
    
    if sizeDeg < sizeCutoffDeg % Small stimulus, RF center varies across days
        goodSessionsAllElectrodes{i} = find(d<=dRange);
    else
        goodSessionsAllElectrodes{i} = 1:numDays;
    end
end
    
osGoodSpikes = []; poGoodSpikes = [];

for i=1:numElectrodes
    disp([i numElectrodes]);
    goodSessionsThisElectrode = goodSessionsAllElectrodes{i};
    numGoodSessions(i) = length(goodSessionsThisElectrode); % Sessions for which stimulus is near the RF
   
    if ~isempty(goodSessionsThisElectrode)
        clear prefOrientationThisElectrode orientationSelectivityThisElectrode frDataThisElectrode 
        for j=1:length(goodSessionsThisElectrode)
            [prefOrientationThisElectrode(j),orientationSelectivityThisElectrode(j),frDataThisElectrode(j,:)] = getOrientationPreferenceFRData(subjectName,expDates{goodSessionsThisElectrode(j)},protocolNames{goodSessionsThisElectrode(j)},folderSourceString,gridType,electrodeList(i),fPos,sPos,frCutoff);
        end

        % GoodSpikes are sessions for which enough spikes were availabe to
        % estimate the OS and PO. For these sessions, OS > -1.
        goodSpkPos = find(orientationSelectivityThisElectrode>-1);
        numGoodSpikeSessions(i) = length(goodSpkPos); % Sessions for which sufficient spikes are available
        osGoodSpikes = [osGoodSpikes orientationSelectivityThisElectrode(goodSpkPos)];
        poGoodSpikes = [poGoodSpikes prefOrientationThisElectrode(goodSpkPos)];

        % GoodOS are sessions for which OS is greater than the osCutoff. We
        % only use these sessions for further analysis.
        goodOSPos = find(orientationSelectivityThisElectrode>osCutoff);
        numGoodOSSessions(i) = length(goodOSPos);
        prefOrientationAllElectrodes{i} = prefOrientationThisElectrode(goodOSPos);
        orientationSelectivityAllElectrodes{i} = orientationSelectivityThisElectrode(goodOSPos);
        
        if ~isempty(goodOSPos)
            frDataAllElectrodes(i,:) = mean(frDataThisElectrode(goodOSPos,:),1);
            [meanPrefOrientation(i),stdPrefOrientation(i),meanOrientationSelectivity(i)] = combineOrientationPreference(prefOrientationAllElectrodes{i},orientationSelectivityAllElectrodes{i});
        else
            meanPrefOrientation(i) = 0;
            stdPrefOrientation(i) = inf;
            meanOrientationSelectivity(i) = 0;
        end
    else
        numGoodSpikeSessions(i) = 0;
        numGoodOSSessions(i) = 0;
        prefOrientationAllElectrodes{i} = [];
        orientationSelectivityAllElectrodes{i} = [];
        
        meanPrefOrientation(i) = 0;
        stdPrefOrientation(i) = inf;
        meanOrientationSelectivity(i) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%% Plot histograms for all data %%%%%%%%%%%%%%%%%%%%%%%
totGoodSessions = sum(numGoodSessions);
disp(['NumElectrodes: ' num2str(numElectrodes) ', sessions: ' num2str(totGoodSessions) ' (' num2str(totGoodSessions/numElectrodes) ' per electrode)']);
totGoodSpikeSessions = sum(numGoodSpikeSessions);
disp(['Sufficient spikes in ' num2str(totGoodSpikeSessions) ' (' num2str(100*totGoodSpikeSessions/totGoodSessions) '%) sessions out of ' num2str(totGoodSessions)]);

osLims = 0.05:0.05:0.95; poLims = 10:10:170;
Nos = hist(osGoodSpikes,osLims);
Npo = hist(poGoodSpikes,poLims);
bar(hOSHistAllElectrodes,osLims,Nos); hold(hOSHistAllElectrodes,'on');
bar(hPOHistAllElectrodes,poLims,Npo); hold(hPOHistAllElectrodes,'on');

totGoodOSSessions = sum(numGoodOSSessions);
disp(['Orientation selectivity higher than cutoff in ' num2str(totGoodOSSessions) ' (' num2str(100*totGoodOSSessions/length(osGoodSpikes)) '%) sessions out of ' num2str(length(osGoodSpikes))]);

badOSPos = find(osGoodSpikes<=osCutoff);
osGoodOS=osGoodSpikes; poGoodOS = poGoodSpikes;
osGoodOS(badOSPos)=[];
poGoodOS(badOSPos)=[];
Nos = hist(osGoodOS,osLims);
Npo = hist(poGoodOS,poLims);
bar(hOSHistAllElectrodes,osLims,Nos,'r'); xlim(hOSHistAllElectrodes,[0 1]);
bar(hPOHistAllElectrodes,poLims,Npo,'r'); xlim(hPOHistAllElectrodes,[0 180]);
xlabel(hOSHistAllElectrodes,'Selectivity');
xlabel(hPOHistAllElectrodes,'Pref Orientation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NumElectrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(hNumElectrodes,electrodeList,numGoodSessions,'go');
hold(hNumElectrodes,'on');
plot(hNumElectrodes,electrodeList,numGoodSpikeSessions,'bv');
plot(hNumElectrodes,electrodeList,numGoodOSSessions,'r*');
xlabel(hNumElectrodes,'Electrode Number');
legend(hNumElectrodes,['All: S=' num2str(totGoodSessions) ',E=' num2str(length(find(numGoodSessions>0)))],['Spikes: S=' num2str(totGoodSpikeSessions) ...
    ',E=' num2str(length(find(numGoodSpikeSessions>0)))],['OS: S=' num2str(totGoodOSSessions) ',E=' num2str(length(find(numGoodOSSessions>0)))],'location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold(hMeanOriAllElectrodes,'on');
hold(hStdPOAllElectrodes,'on');
        
for i=1:numElectrodes
    if numGoodOSSessions(i)>0
        if stdPrefOrientation(i)<stdPOCutoff
            colorMarkerAllStr = 'co';
            colorMarkerMeanStr = 'kV';
        else
            colorMarkerAllStr = 'go';
            colorMarkerMeanStr = 'rV';
        end
        plot(hMeanOriAllElectrodes,electrodeList(i)+zeros(1,numGoodOSSessions(i)),prefOrientationAllElectrodes{i},colorMarkerAllStr(1));
        plot(hMeanOriAllElectrodes,electrodeList(i),prefOrientationAllElectrodes{i},colorMarkerAllStr);
        plot(hMeanOriAllElectrodes,electrodeList(i),meanPrefOrientation(i),colorMarkerMeanStr);
        
        plot(hStdPOAllElectrodes,electrodeList(i),stdPrefOrientation(i),colorMarkerMeanStr);
    end
end

ylabel(hMeanOriAllElectrodes,'Preferred Orientation (Deg)');
ylabel(hStdPOAllElectrodes,'Circ Std Dev (between 0 and 1)');
xlabel(hStdPOAllElectrodes,'Electrode Number');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% histograms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finalGoodPos = find(stdPrefOrientation<stdPOCutoff);

finalElectrodeList = electrodeList(finalGoodPos);
finalOrientationPref = meanPrefOrientation(finalGoodPos);
finalOrientationSelectivity = meanOrientationSelectivity(finalGoodPos);
finalFiringRates = frDataAllElectrodes(finalGoodPos,:);

numFinalElectrodes = length(finalElectrodeList);
disp([num2str(numFinalElectrodes) ' have reliable Preferred orientation']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Histograms %%%%%%%%%%%%%%%%%%%%%%%%%%%
stdLims = osLims;
Nstd = hist(stdPrefOrientation,stdLims);
bar(hStdPOHistAllElectrodes,stdLims,Nstd,'r');
hold(hStdPOHistAllElectrodes,'on');

NstdFinal = hist(stdPrefOrientation(finalGoodPos),stdLims);
bar(hStdPOHistAllElectrodes,stdLims,NstdFinal,'k'); xlim(hStdPOHistAllElectrodes,[0 1]);

NpoFinal = hist(meanPrefOrientation(finalGoodPos),poLims);
bar(hPOHistAllElectrodes,poLims,NpoFinal,'k'); xlim(hPOHistAllElectrodes,[0 180]);

%%%%%%%%%%%%%%%%%% Display Orientations in a color map %%%%%%%%%%%%%%%%%%%%
colorMapOrientations = hsv(180);

% Show ColorMap
% hold(hColorMapPlot,'on');
% for i=1:3:180
%     for j=0:0.1:1
%         plot(hColorMapPlot,j*cosd(i),j*sind(i),'marker','o','color',j*colorMapOrientations(i,:),'markerFaceColor',j*colorMapOrientations(i,:));
%     end
% end
% axis(hColorMapPlot,[-1 1 -1 1]);

hold(hPORFPlotAllElectrodes,'on');
hold(hOSRFPlotAllElectrodes,'on');

for i=1:numElectrodes
    a = rfData.rfStats(electrodeList(i)).meanAzi;
    e = rfData.rfStats(electrodeList(i)).meanEle;
    
    if (stdPrefOrientation(i)>stdPOCutoff)
        plot(hPORFPlotAllElectrodes,a,e,'marker','o','color',[0 0 0]);
        plot(hOSRFPlotAllElectrodes,a,e,'marker','o','color',[0 0 0]);
    else
        colorVal = colorMapOrientations(ceil(meanPrefOrientation(i)),:);
        plot(hPORFPlotAllElectrodes,a,e,'marker','o','markersize',12,'color',colorVal,'markerFaceColor',colorVal);
        plot(hOSRFPlotAllElectrodes,a,e,'marker','o','markersize',12,'color',[1 1 1]*meanOrientationSelectivity(i),'markerFaceColor',[1 1 1]*meanOrientationSelectivity(i));
        
        textOriStr = [num2str(round(meanPrefOrientation(i))) '\circ'];
        text(a,e,textOriStr,'parent',hPORFPlotAllElectrodes,'HorizontalAlignment','center');
        textOSStr = num2str(round(meanOrientationSelectivity(i),2));
        text(a,e,textOSStr,'parent',hOSRFPlotAllElectrodes,'HorizontalAlignment','center');
    end
end

axis(hPORFPlotAllElectrodes,gridLims);
axis(hOSRFPlotAllElectrodes,gridLims);

soValsUnique = 0:0.2:1;
for i=1:length(oValsUnique) 
    text(gridLims(1)+0.1,gridLims(4)-0.15*i,[num2str(oValsUnique(i)) '\circ'],'color',colorMapOrientations(oValsUnique(i)+1,:),'parent',hPORFPlotAllElectrodes,'fontsize',14,'HorizontalAlignment','center');
    text(gridLims(1)+0.1,gridLims(4)-0.15*i,num2str(soValsUnique(i)),'color',[1 1 1]*soValsUnique(i),'parent',hOSRFPlotAllElectrodes,'fontsize',14,'HorizontalAlignment','center');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefOrientation,stdPrefOrientation,meanOrientationSelectivity] = combineOrientationPreference(prefOrientationThisElectrode,orientationSelectivityThisElectrode)

% Orientations vary between 0 to 180, with pref of 180 degrees very close
% to 0. To compute the circular mean, we multiply orientations by 2, take
% the mean, and then divide by 2.

poInRad = (2*pi)*(prefOrientationThisElectrode(:)/180); % These orientations vary between 0 and 2pi.
mpoInRad = circ_mean(poInRad); % This is between -pi and pi
tmp = find(mpoInRad<0);
mpoInRad(tmp) = mpoInRad(tmp) + 2*pi; % This is between 0 and 2pi
prefOrientation = 180*mpoInRad/(2*pi); % These are between -180 and 180.

stdPrefOrientation = circ_std(poInRad(:))/sqrt(2); % Scaled so that it varies between 0 and 1.
meanOrientationSelectivity = mean(orientationSelectivityThisElectrode);
end
function displaySavedRFcenters(hGridPlot,plotHandle,subjectName,gridType,electrodeList,axisLims)

rfData = load([subjectName gridType 'RFData.mat']);
rfStats = rfData.rfStats;

startRow = 1;
[~,electrodeColorNames,electrodeArray] = showElectrodeLocationsInColor([],hGridPlot,1,[],startRow,gridType,subjectName);

% Only show the electrodes in the list
showElectrodeLocations([],setdiff(1:96,electrodeList),'w',hGridPlot,1);

for i=1:length(electrodeList)
 [row,col] = find(electrodeList(i) == electrodeArray);
 colorList{i} = electrodeColorNames{row,col}; %#ok<*AGROW>
end

set(plotHandle,'NextPlot','add');
for i=1:length(electrodeList)
    rfParams = rfStats(electrodeList(i)).params;
    plot(plotHandle,rfParams(1),rfParams(2),'color',colorList{i},'Marker','O','markersize',10);
    textENum = num2str(electrodeList(i));
    text(rfParams(1),rfParams(2),textENum,'parent',plotHandle,'HorizontalAlignment','center');
end

set(plotHandle,'NextPlot','replace');
axis(plotHandle,axisLims);
end
function [prefOrientation,orientationSelectivity,fr,oValsUnique] = getOrientationPreferenceFRData(subjectName,expDate,protocolName,folderSourceString,gridType,electrodeNum,fPos,sigmaPos,frCutoff)

timeRangeForFRComputation = [0 0.2];
delT = diff(timeRangeForFRComputation);
stimPosNum=3;

[~,~,~,~,oValsUnique,~,~,parameterCombinations] = loadParameterCombinations(subjectName,expDate,protocolName,folderSourceString,gridType);

% Get bad trials
load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'));
[allStimPos{1},allStimPos{2}]= getStimPos(subjectName,expDate,protocolName,folderSourceString,gridType);
allStimPos{3} = unique([allStimPos{1} allStimPos{2}]);

% Spike Data
spikeDataFile = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','Spikes',['elec' num2str(electrodeNum) '_SID0.mat']);
if exist(spikeDataFile,'file')
    load(spikeDataFile);
else
    spikeData=[];
end

for i=1:length(oValsUnique)
    goodPos = intersect(parameterCombinations{1,1,sigmaPos,fPos,i},allStimPos{stimPosNum});
    goodPos = setdiff(goodPos,badTrials);
    if isempty(spikeData)
        fr(i) = 0;
    else
        fr(i) = mean(getSpikeCounts(spikeData(goodPos),timeRangeForFRComputation))/delT;
    end
end

if max(fr)<= frCutoff % Max firing rate should be at least frCutoff spikes/s
    prefOrientation = 0;
    orientationSelectivity = -1;
else
    num=0;
    den=0;
    
    for i=1:length(oValsUnique)
        num = num+fr(i)*sind(2*oValsUnique(i));
        den = den+fr(i)*cosd(2*oValsUnique(i));
    end
    
    prefOrientation = 90*atan2(num,den)/pi;
    orientationSelectivity = abs(den+1i*num)/sum(fr);
end

tmp = find(prefOrientation<0);
prefOrientation(tmp) = prefOrientation(tmp)+180;

end
function [allStimPosEqualToOne,allStimPosGreaterThanOne]= getStimPos(subjectName,expDate,protocolName,folderSourceString,gridType)

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');

load(fullfile(folderExtract,'goodStimNums.mat'));
load(fullfile(folderExtract,'stimResults.mat'));
goodStimPos = stimResults.stimPosition(goodStimNums);

if exist([folderExtract 'validStimAfterTarget.mat'],'file')
    load([folderExtract 'validStimAfterTarget.mat']);
    goodStimPos(validStimuliAfterTarget)=-1;  % These will be not be included in either stimPos==1 or stimPos>1
end

allStimPosEqualToOne = find(goodStimPos==1);
allStimPosGreaterThanOne = find(goodStimPos>1);
end
function [aValsUnique,eValsUnique,srValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique,parameterCombinations,srTag] = loadParameterCombinations(subjectName,expDate,protocolName,folderSourceString,gridType) %#ok<*STOUT>
load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat'));

if exist('rValsUnique','var')
    srValsUnique=rValsUnique;
    srTag = 'r';
elseif exist('sValsUnique','var')
    srValsUnique=sValsUnique;
    srTag = 's';
end

if ~exist('cValsUnique','var')
    cValsUnique=100;
end
if ~exist('tValsUnique','var')
    tValsUnique=0;
end
end