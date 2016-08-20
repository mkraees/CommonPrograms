% Modified from showEyeDataRafikiDiffConditionsGRFLongTrials

function displayEyeDataGRFLongTrials(subjectName,expDates,protocolNames,folderSourceString,gridType,protocolType)

if ~iscell(expDates)
    [expDates,protocolNames] = convertToCell(expDates,protocolNames);
end

%%%% display properties
fontSizeSmall=10;
fontSizeMedium=12;
fontSizeLarge=14;
backgroundColor = 'w';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots 

hEyePlots        = getPlotHandles(1,2,[0.05 0.60 0.4 0.3],0.05,0);
hEyeSigPlots     = getPlotHandles(1,2,[0.05 0.50 0.4 0.1],0.05,0);

hEyeMeanPlots    = getPlotHandles(2,1,[0.5 0.50 0.10 0.4],0.05,0);
hEyeMeanSigPlots = getPlotHandles(2,1,[0.6 0.50 0.10 0.4],0.05,0);
hEyePlotVsStimPos = getPlotHandles(2,1,[0.75 0.50 0.20 0.4],0.05,0);

hEyeMeanPlotsSingleDay = getPlotHandles(2,1,[0.75 0.05 0.2 0.4],0.05,0);

hMSPlot         = subplot('position',[0.275 0.05 0.175 0.4]);
%hMSSigPlot      = subplot('position',[0.05 0.05 0.175 0.1]);
hSpeedHistogram  = subplot('position',[0.05 0.05 0.175 0.4]);

hMSMeanPlot      = subplot('position',[0.5 0.05 0.10 0.4]);
hMSMeanSigPlot   = subplot('position',[0.6 0.05 0.10 0.4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls

hChoicePanel = uipanel('Title','Day','Unit','Normalized','fontSize', fontSizeLarge,'Position',[0 0.925 0.15 0.075]);

% Make protocolString
protocolString = getProtocolString(subjectName,expDates,protocolNames,folderSourceString,gridType);
uicontrol('Parent',hChoicePanel,'Unit','Normalized', 'Position',[0 0.5 0.5 0.5], ...
    'Style','text','String','Day','FontSize',fontSizeSmall);
hDay = uicontrol('Parent',hChoicePanel,'Unit','Normalized', 'Position',[0.5 0.5 0.5 0.5], ...
    'BackgroundColor', backgroundColor, ...
    'Style','popup','String',protocolString,'FontSize',fontSizeSmall);

% Speed cut off ranges
uicontrol('Parent',hChoicePanel,'Unit','Normalized', 'Position',[0 0 0.5 0.5], ...
    'Style','text','String','Cutoff(Deg/s):','FontSize',fontSizeSmall);
hCutoff = uicontrol('Parent',hChoicePanel,'Unit','Normalized', 'Position',[0.5 0 0.5 0.5], ...
    'BackgroundColor', backgroundColor, ...
    'Style','edit','String','15','FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot, Rescale, Hold and cla
hPlottingPanel = uipanel('Title','Plot options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.15 0.925 0.15 0.075]);

uicontrol('Parent',hPlottingPanel,'Unit','Normalized', 'Position',[0 0.5 0.5 0.5], ...
    'Style','pushbutton','String','Plot','FontSize',fontSizeMedium, ...
    'Callback',{@plot_Callback});
% uicontrol('Parent',hPlottingPanel,'Unit','Normalized', 'Position',[0.5 0.5 0.5 0.5], ...
%     'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium, ...
%     'Callback',{@rescale_Callback});
% uicontrol('Parent',hPlottingPanel,'Unit','Normalized','Position',[0 0 0.5 0.5], ...
%     'Style','pushbutton','String','plot size','FontSize',fontSizeMedium, ...
%     'Callback',{@replotSizePlot_Callback});
uicontrol('Parent',hPlottingPanel,'Unit','Normalized', 'Position',[0.5 0.5 0.5 0.5], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

uicontrol('Parent',hPlottingPanel,'Unit','Normalized', 'Position',[0 0 0.5 0.5], ...
    'Style','text','String','Test Method','FontSize',fontSizeSmall);

[testMethodString,testMethods] = getTestMethods;
hTestMethod = uicontrol('Parent',hPlottingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position', [0.5 0 0.5 0.5], ...
    'Style','popup','String',testMethodString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the parameters. They should be the same for all the days
folderExtract = fullfile(folderSourceString,'data',subjectName,gridType,expDates{1},protocolNames{1},'extractedData');
[~,~,~,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);

parameterTextWidth = 0.2; parameterWidth = 0.1;
hParameterPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.3 0.925 0.4 0.075]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuth
% %azimuthString = getStringFromValues(aValsUnique,1);
% uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
%     'Position',[0 0.5 parameterTextWidth 0.5],...
%     'Style','text','String','Azimuth (Deg)','FontSize',fontSizeSmall);
% hAzimuth = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
%     'Position', [parameterTextWidth 0.5 parameterWidth 0.5], ...
%     'Style','text','String','','FontSize',fontSizeSmall);
% 
% % Elevation
% %elevationString = getStringFromValues(eValsUnique,1);
% uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
%     'Position',[0 0 parameterTextWidth 0.5], ...
%     'Style','text','String','Elevation (Deg)','FontSize',fontSizeSmall);
% hElevation = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
%     'Position', [parameterTextWidth 0 parameterWidth 0.5], ...
%     'Style','text','String','','FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 0.5 parameterTextWidth 0.5], ...
    'Style','text','String','Sigma (Deg)','FontSize',fontSizeSmall);

sigmaString = getStringFromValues(sValsUnique,1);
hSigma = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position', [parameterTextWidth 0.5 parameterWidth 0.5], ...
    'Style','popup','String',sigmaString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial Frequency
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 0 parameterTextWidth 0.5], ...
    'Style','text','String','SF (CPD)','FontSize',fontSizeSmall);

spatialFreqString = getStringFromValues(fValsUnique,1);
hSpatialFreq = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position', [parameterTextWidth 0 parameterWidth 0.5], ...
    'Style','popup','String',spatialFreqString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orientation
orientationString = getStringFromValues(oValsUnique,1);
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[1*(parameterTextWidth+parameterWidth) 0.5 parameterTextWidth 0.5], ...
    'Style','text','String','Orientation (Deg)','FontSize',fontSizeSmall);
hOrientation = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [1*(parameterTextWidth+parameterWidth)+parameterTextWidth 0.5 parameterWidth 0.5], ...
    'Style','popup','String',orientationString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contrast
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[1*(parameterTextWidth+parameterWidth) 0 parameterTextWidth 0.5], ...
    'Style','text','String','Contrast (%)','FontSize',fontSizeSmall);

contrastString = getStringFromValues(cValsUnique,1);
hContrast = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [1*(parameterTextWidth+parameterWidth)+parameterTextWidth 0 parameterWidth 0.5], ...
    'Style','popup','String',contrastString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal Frequency
temporalFreqString = getStringFromValues(tValsUnique,1);
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[2*(parameterTextWidth+parameterWidth) 0 parameterTextWidth 0.5], ...
    'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeSmall);
hTemporalFreq = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [2*(parameterTextWidth+parameterWidth)+parameterTextWidth 0 parameterWidth 0.5], ...
    'Style','popup','String',temporalFreqString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overwrite
if strcmp(protocolType,'SF') % Compare SF
    set(hSpatialFreq,'Style','text','String','tested','FontSize',fontSizeSmall);

elseif strcmp(protocolType,'TF')
    set(hTemporalFreq,'Style','text','String','tested','FontSize',fontSizeSmall);

elseif strcmp(protocolType,'CON')
    set(hContrast,'Style','text','String','tested','FontSize',fontSizeSmall);
    
else                        % compare Sigma
    set(hSigma,'Style','text','String','tested','FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Timing controls
hTimingPanel = uipanel('Title','Timing','Unit','Normalized','fontSize', fontSizeLarge, 'Position',[0.7 0.925 0.3 0.075]);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 0.5 0.25 0.5], ...
    'Style','text','String','Time Period','FontSize',fontSizeSmall);

tRange = [0.2 0.4];
hTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.25 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(tRange(1)),'FontSize',fontSizeSmall);

hTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(tRange(2)),'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    function cla_Callback(~,~)
        
        claGivenPlotHandle(hEyePlots); claGivenPlotHandle(hEyeSigPlots);
        claGivenPlotHandle(hEyePlotVsStimPos);
        claGivenPlotHandle(hEyeMeanPlots); claGivenPlotHandle(hEyeMeanSigPlots);
        claGivenPlotHandle(hEyeMeanPlotsSingleDay);
        
        claGivenPlotHandle(hMSMeanPlot); claGivenPlotHandle(hMSMeanSigPlot);
        claGivenPlotHandle(hMSPlot); claGivenPlotHandle(hSpeedHistogram);
    end
    function plot_Callback(~,~)
        
        claGivenPlotHandle(hEyePlots); claGivenPlotHandle(hEyeSigPlots);
        claGivenPlotHandle(hEyePlotVsStimPos);
        claGivenPlotHandle(hEyeMeanPlots); claGivenPlotHandle(hEyeMeanSigPlots);
        
        claGivenPlotHandle(hMSMeanPlot); claGivenPlotHandle(hMSMeanSigPlot);
        claGivenPlotHandle(hMSPlot); claGivenPlotHandle(hSpeedHistogram);
        %claGivenPlotHandle(hEyeMeanPlotsSingleDay);
        
        dayNum = get(hDay,'val');
        cutoff = str2double(get(hCutoff,'string'));
        useThisTestMethod = testMethods{get(hTestMethod,'val')};
        
        useThisTimeRange=[str2double(get(hTMin,'string')) str2double(get(hTMax,'string'))];
        
        %Initialze
        clear expDate protocolName
        if dayNum > length(expDates)
            expDate = expDates;
            protocolName = protocolNames;
        else
            expDate = expDates{dayNum};
            protocolName = protocolNames{dayNum};
        end
        
        % a, e and o are the same for all protocolListNumbers
        a=1;e=1;o=get(hOrientation,'val');
        
        if strcmp(protocolType,'SF') % SF will be compared
            numberOfCategories = length(fValsUnique);
            s = get(hSigma,'val');
            f = numberOfCategories+1;  % Not required, computed for completeness
            c = get(hContrast,'val'); 
            t = get(hTemporalFreq,'val');
            
        elseif strcmp(protocolType,'TF') % TF will be compared
            numberOfCategories = length(tValsUnique);
            s = get(hSigma,'val');
            f = get(hSpatialFreq,'val');
            c = get(hContrast,'val');  
            t = numberOfCategories+1; % Not required, computed for completeness
        
        elseif strcmp(protocolType,'CON') % Contrast will be compared
            numberOfCategories = length(cValsUnique);
            s = get(hSigma,'val');
            f = get(hSpatialFreq,'val');
            c = numberOfCategories+1; % Not required, computed for completeness 
            t = get(hTemporalFreq,'val');
            
        elseif strcmp(protocolType,'ORI') % Contrast will be compared
            numberOfCategories = length(oValsUnique);
            s = get(hSigma,'val');
            f = get(hSpatialFreq,'val');
            c = get(hContrast,'val');
            t = get(hTemporalFreq,'val');
            o = numberOfCategories+1; % Not required, computed for completeness
            
        else % Sigma will be compared
            numberOfCategories = length(sValsUnique);
            s = numberOfCategories+1; % Not required, computed for completeness
            f = get(hSpatialFreq,'val');
            c = get(hContrast,'val'); 
            t = get(hTemporalFreq,'val');
        end
        
        clear eyeX eyeY
        colorNames = jet(numberOfCategories);
        [eyeX,eyeY,xs,~,~,comparisonDataValues,MSData,allEyeSpeeds,numMSInRange] = getSortedEyeData(subjectName,expDate,protocolName,folderSourceString,gridType,protocolType,a,e,s,f,o,c,t,numberOfCategories,cutoff,useThisTimeRange);

        % Updat the azimuth and elevation fields
        %set(hAzimuth,'String',aziStr); set(hElevation,'String',eleStr);
        
        % Display Eye Data into categories
        compareAndDisplayData(eyeX,xs,hEyePlots(1),hEyeSigPlots(1),useThisTestMethod,colorNames);%,colorX,colorSig,displayPlot,showStdErr);
        compareAndDisplayData(eyeY,xs,hEyePlots(2),hEyeSigPlots(2),useThisTestMethod,colorNames);
        
        % Put labels
        ylabel(hEyePlots(1),'x eyedata (deg)'); ylabel(hEyePlots(2),'y eyedata (deg)');
        xlabel(hEyeSigPlots(1),'time (s)'); ylabel(hEyeSigPlots(1),'p-value');       
        xlabel(hEyeSigPlots(2),'time (s)'); ylabel(hEyeSigPlots(2),'p-value'); 
        
        %put legend
        numComparisonDataValues = length(comparisonDataValues);
        legendStr1 = cell(1,numComparisonDataValues);
        legendStr2 = cell(1,numComparisonDataValues);
        for legendPos=1:numComparisonDataValues
            legendStr1{legendPos} = num2str(comparisonDataValues(legendPos));
            legendStr2{legendPos} =  ['n=' num2str(size(eyeX{legendPos},1))];
        end
        legend(hEyePlots(1),legendStr1,'Location','NorthOutside');
        legend(hEyePlots(2),legendStr2,'Location','NorthOutside');     
        
        % Significance analysis on Mean Data
        compareAndDisplayMeanData(eyeX,xs,hEyeMeanPlots(1),hEyeMeanSigPlots(1),useThisTestMethod,useThisTimeRange,colorNames);
        compareAndDisplayMeanData(eyeY,xs,hEyeMeanPlots(2),hEyeMeanSigPlots(2),useThisTestMethod,useThisTimeRange,colorNames);
        ylabel(hEyeMeanPlots(1),'x position (deg)'); ylabel(hEyeMeanPlots(2),'y position (deg)');
        
        % Plot eye position versus stimulus Pos
        plotEyePosititionVsStimPos(subjectName,expDate,protocolName,folderSourceString,gridType,hEyePlotVsStimPos(1),hEyePlotVsStimPos(2));

        % Show Microsaccades      
        % Histogram
        eyeSpeedCenters=0:100;
        histAllSpeeds = hist(allEyeSpeeds,eyeSpeedCenters);
        plot(hSpeedHistogram,eyeSpeedCenters,log10(histAllSpeeds),'k');
        hold(hSpeedHistogram,'on');
        MSPos = find(eyeSpeedCenters>cutoff);
        plot(hSpeedHistogram,eyeSpeedCenters(MSPos),log10(histAllSpeeds(MSPos)),'r');
        hold(hSpeedHistogram,'off');
        ylabel(hSpeedHistogram,'log(#)');
        xlabel(hSpeedHistogram,'speed (deg/s)');
        
 
        % # microsaccades /sec
        for ii=1:numberOfCategories
            [H,timeValsMS] = getPSTH(MSData{ii},20,[xs(1) xs(end)]);
            plot(hMSPlot,timeValsMS,H,'color',colorNames(ii,:));
            hold(hMSPlot,'on');
        end
        hold(hMSPlot,'off');
        axis(hMSPlot,'tight');
        xlabel(hMSPlot,'time (s)'); ylabel(hMSPlot,'MicroSaccades/s');
        
        % Significance test
        compareAndDisplayMeanData2(numMSInRange,hMSMeanPlot,hMSMeanSigPlot,useThisTestMethod,colorNames);
        
         % if expDates is a cell array (population data, report the p-values
        % of individual days as well.
        
        if iscell(expDate)
            K=length(expDate);
            for ii=1:K
                clear eyeX eyeY xs aziStr eleStr comparisonDataValues
                disp(['Performing significance analysis on ' num2str(ii) ': ' expDate{ii} protocolName{ii}]);
                [eyeX,eyeY,xs] = getSortedEyeData(subjectName,expDate{ii},protocolName{ii},folderSourceString,gridType,protocolType,a,e,s,f,o,c,t,numberOfCategories,cutoff,useThisTimeRange);
                pXAllDays(ii)=compareAndDisplayMeanData(eyeX,xs,[],[],useThisTestMethod,useThisTimeRange,[],0); %#ok<*AGROW>
                pYAllDays(ii)=compareAndDisplayMeanData(eyeY,xs,[],[],useThisTestMethod,useThisTimeRange,[],0);
                
                for iii=1:length(eyeX)
                    numberOfStimuliPerCondition(ii,iii) = size(eyeX{iii},1);
                end
            end
           
            disp(['Mean number of stimuli per condition: ' num2str(mean(mean(numberOfStimuliPerCondition)))]);
            disp(['Max number of stimuli per condition: ' num2str(max(max(numberOfStimuliPerCondition)))]);
            disp(['Min number of stimuli per condition: ' num2str(min(min(numberOfStimuliPerCondition)))]);
            
            plotSignificanceData(pXAllDays,1:K,hEyeMeanPlotsSingleDay(1),'k',0.05);
            plotSignificanceData(pYAllDays,1:K,hEyeMeanPlotsSingleDay(2),'k',0.05);
            
            ylabel(hEyeMeanPlotsSingleDay(1),'p-value X data');
            ylabel(hEyeMeanPlotsSingleDay(2),'p-value Y data');
            xlabel(hEyeMeanPlotsSingleDay(2),'day number');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function protocolString = getProtocolString(subjectName,expDates,protocolNames,folderSourceString,gridType)
protocolString='';
nTotal=0;
for i=1:length(expDates)
    expDate = expDates{i};
    protocolName = protocolNames{i};
    
    % Find out the number of stimuli
    folderExtract = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData');
    N = getNumStimuli(folderExtract);
    protocolString = cat(2,protocolString,[num2str(i) ' ' expDate protocolName '(N=' num2str(N) ')|']);
    nTotal=nTotal+N;
end

protocolString = [protocolString 'all Days (N=' num2str(nTotal) ')'];
end
function outString = getStringFromValues(valsUnique,decimationFactor)

if length(valsUnique)==1
    outString = convertNumToStr(valsUnique(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique)
        outString = cat(2,outString,[convertNumToStr(valsUnique(i),decimationFactor) '|']);
    end
    outString = [outString 'all'];
end

    function str = convertNumToStr(num,f)
        if num > 16384
            num=num-32768;
        end
        str = num2str(num/f);
    end
end
function claGivenPlotHandle(plotHandles)
[numRows,numCols] = size(plotHandles);
for ii=1:numRows
    for jj=1:numCols
        cla(plotHandles(ii,jj));
    end
end
end
function numStimuli = getNumStimuli(folderExtract)
clear goodStimNums
load(fullfile(folderExtract,'goodStimNums.mat'));

if exist([folderExtract 'validStimAfterTarget.mat'],'file')
    load(fullfile(folderExtract,'validStimAfterTarget.mat'));
    %disp(['Removing ' num2str(length(validStimuliAfterTarget)) ' stimuli from goodStimNums']);
    goodStimNums(validStimuliAfterTarget)=[];
end

clear stimResults
load(fullfile(folderExtract,'stimResults.mat'));

goodStimPos = stimResults.stimPosition(goodStimNums);

% all stimPostinios greater than 1
% numStimuli = length(find(goodStimPos>1));
numStimuli = length(find(goodStimPos>0)); % All stimuli, including the first one
end
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract) %#ok<*STOUT>

load(fullfile(folderExtract,'parameterCombinations.mat'));

if      ~exist('cValsUnique','var'),  cValsUnique=[];                   end
if      ~exist('tValsUnique','var'),  tValsUnique=[];                   end
if      ~exist('sValsUnique','var'),  sValsUnique=rValsUnique/3;        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eyeX,eyeY,timeValsEyePos,aziStr,eleStr,comparisonDataValues,MSData,allEyeSpeeds,numMSInRange] = getSortedEyeData(subjectName,expDates,protocolNames,folderSourceString,gridType,protocolType,a,e,s,f,o,c,t,numberOfCategories,cutoff,timeRange)

FsEye=200;

if ~iscell(expDates)
    [expDates,protocolNames] = convertToCell(expDates,protocolNames);
end

% Initialization
eyeX=cell(1,numberOfCategories);
eyeY=cell(1,numberOfCategories);
MSData=cell(1,numberOfCategories);
numMSInRange=cell(1,numberOfCategories);

allEyeSpeeds=[];

for i=1:length(expDates)
    expDate = expDates{i};
    protocolName=protocolNames{i};
    
    folderExtract = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData');
    folderSegment = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData');
    
    % get timevals for eye position. 
    if i==1
        load(fullfile(folderExtract,'EyeData.mat')); % returns the variable 'timeVals'
        timeValsEyePos = (eyeRangeMS(1):1000/FsEye:eyeRangeMS(2)-1000/FsEye)/1000;
    end
    
    % The stimuli for which eyeData has been stored (stimPos>1)
    clear goodStimNums
    load(fullfile(folderExtract,'goodStimNums.mat'));
    
    clear stimResults
    load(fullfile(folderExtract,'stimResults.mat'));

    clear goodStimPos useTheseStims
    goodStimPos = stimResults.stimPosition(goodStimNums);

    if exist([folderExtract 'validStimAfterTarget.mat'],'file')
        load([folderExtract 'validStimAfterTarget.mat']);
        goodStimPos(validStimuliAfterTarget)=-1;
    end

    % all stimPostions greater than 1
    % allUsefulStims = find(goodStimPos>1);
    allUsefulStims = find(goodStimPos>0);
    
    % eyeData
    load(fullfile(folderSegment,'eyeData','eyeDataDeg'));
    
    % eyeSpeed
    eyeDataDegX=cell2mat(eyeDataDegX)';
    eyeDataDegY=cell2mat(eyeDataDegY)';
    lengthEyeSignal = size(eyeDataDegX,2);
    
    for j=1:size(eyeDataDegX,1)
        eyeSpeedX(j,:) = [eyeDataDegX(j,2:lengthEyeSignal)-eyeDataDegX(j,1:lengthEyeSignal-1) 0];
        eyeSpeedY(j,:) = [eyeDataDegY(j,2:lengthEyeSignal)-eyeDataDegY(j,1:lengthEyeSignal-1) 0];
    end
    
    eyeSpeedMag = FsEye*sqrt(eyeSpeedX.^2+eyeSpeedY.^2);
    
    if size(eyeDataDegX,1) ~= length(allUsefulStims) %#ok<*NODEF>
        error('eyeData and allUsefulStims have unequal lengths');
    else
        inverseMap=zeros(1,length(goodStimPos));
        inverseMap(allUsefulStims)=1:length(allUsefulStims);
    end
    
    % parameterCombinations
    [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);

    if strcmp(protocolType,'SF')
        
        for j=1:numberOfCategories
            eyeIndices = inverseMap(intersect(parameterCombinations{a,e,s,j,o,c,t},allUsefulStims));
            eyeX{j} = cat(1,eyeX{j},eyeDataDegX(eyeIndices,:));
            eyeY{j} = cat(1,eyeY{j},eyeDataDegY(eyeIndices,:));
            
            clear MStmp numMStmp
            [MStmp,numMStmp] = findMicroSaccades(eyeSpeedMag(eyeIndices,:),cutoff,timeValsEyePos,timeRange);
            MSData{j} = [MSData{j} MStmp];
            numMSInRange{j}  = [numMSInRange{j} numMStmp];
            
            clear tmpEyeSpeeds
            tmpEyeSpeeds = eyeSpeedMag(eyeIndices,:);
            allEyeSpeeds=cat(1,allEyeSpeeds,tmpEyeSpeeds(:));
        end
        
        comparisonDataValues = fValsUnique;
        
    elseif strcmp(protocolType,'TF')
        
        for j=1:numberOfCategories
            eyeIndices = inverseMap(intersect(parameterCombinations{a,e,s,f,o,c,j},allUsefulStims));
            eyeX{j} = cat(1,eyeX{j},eyeDataDegX(eyeIndices,:));
            eyeY{j} = cat(1,eyeY{j},eyeDataDegY(eyeIndices,:));
            
            clear MStmp numMStmp
            [MStmp,numMStmp] = findMicroSaccades(eyeSpeedMag(eyeIndices,:),cutoff,timeValsEyePos,timeRange);
            MSData{j} = [MSData{j} MStmp];
            numMSInRange{j}  = [numMSInRange{j} numMStmp];
            
            clear tmpEyeSpeeds
            tmpEyeSpeeds = eyeSpeedMag(eyeIndices,:);
            allEyeSpeeds=cat(1,allEyeSpeeds,tmpEyeSpeeds(:));
        end
        
        comparisonDataValues = tValsUnique;

    elseif strcmp(protocolType,'CON')
        
        for j=1:numberOfCategories
            eyeIndices = inverseMap(intersect(parameterCombinations{a,e,s,f,o,j,t},allUsefulStims));
            eyeX{j} = cat(1,eyeX{j},eyeDataDegX(eyeIndices,:));
            eyeY{j} = cat(1,eyeY{j},eyeDataDegY(eyeIndices,:));
            
            clear MStmp numMStmp
            [MStmp,numMStmp] = findMicroSaccades(eyeSpeedMag(eyeIndices,:),cutoff,timeValsEyePos,timeRange);
            MSData{j} = [MSData{j} MStmp];
            numMSInRange{j}  = [numMSInRange{j} numMStmp];
            
            clear tmpEyeSpeeds
            tmpEyeSpeeds = eyeSpeedMag(eyeIndices,:);
            allEyeSpeeds=cat(1,allEyeSpeeds,tmpEyeSpeeds(:));
        end
        
        comparisonDataValues = cValsUnique;
        
    elseif strcmp(protocolType,'ORI')
        
        for j=1:numberOfCategories
            eyeIndices = inverseMap(intersect(parameterCombinations{a,e,s,f,o,j,t},allUsefulStims));
            eyeX{j} = cat(1,eyeX{j},eyeDataDegX(eyeIndices,:));
            eyeY{j} = cat(1,eyeY{j},eyeDataDegY(eyeIndices,:));
            
            clear MStmp numMStmp
            [MStmp,numMStmp] = findMicroSaccades(eyeSpeedMag(eyeIndices,:),cutoff,timeValsEyePos,timeRange);
            MSData{j} = [MSData{j} MStmp];
            numMSInRange{j}  = [numMSInRange{j} numMStmp];
            
            clear tmpEyeSpeeds
            tmpEyeSpeeds = eyeSpeedMag(eyeIndices,:);
            allEyeSpeeds=cat(1,allEyeSpeeds,tmpEyeSpeeds(:));
        end
        
        comparisonDataValues = oValsUnique;
    else
        
        for j=1:numberOfCategories
            eyeIndices = inverseMap(intersect(parameterCombinations{a,e,j,f,o,c,t},allUsefulStims));
            eyeX{j} = cat(1,eyeX{j},eyeDataDegX(eyeIndices,:));
            eyeY{j} = cat(1,eyeY{j},eyeDataDegY(eyeIndices,:));
            
            clear MStmp numMStmp
            [MStmp,numMStmp] = findMicroSaccades(eyeSpeedMag(eyeIndices,:),cutoff,timeValsEyePos,timeRange);
            MSData{j} = [MSData{j} MStmp];
            numMSInRange{j}  = [numMSInRange{j} numMStmp];
            
            clear tmpEyeSpeeds
            tmpEyeSpeeds = eyeSpeedMag(eyeIndices,:);
            allEyeSpeeds=cat(1,allEyeSpeeds,tmpEyeSpeeds(:));
        end
        
        comparisonDataValues = sValsUnique;
    end
end

if i==1
    aziStr = num2str(aValsUnique(1));
    eleStr = num2str(eValsUnique(1));
else
    aziStr = 'many';
    eleStr = 'many';
end
end
function [testMethodString,testMethods] = getTestMethods

testMethods{1} = 'default';
testMethods{2} = 'anova';
testMethods{3} = 'kruskalWallis';
testMethods{4} = 'ttest';

testMethodString='';
for i=1:length(testMethods)
    testMethodString = cat(2,testMethodString,[testMethods{i} '|']);
end
testMethodString = removeIfPresent(testMethodString,'|');
end
function [expDates,protocolNames] = convertToCell(expDate,protocolName)
expDates{1} = expDate;
protocolNames{1} = protocolName;
end
