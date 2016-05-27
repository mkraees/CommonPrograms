% Display All Channels
function displayAllChannelsGRF(subjectName,expDate,protocolName,folderSourceString,gridType,gridLayout)

if ~exist('folderSourceString','var');   folderSourceString='E:';       end
if ~exist('gridType','var');             gridType='EEG';                end
if ~exist('gridLayout','var');          gridLayout=0;                   end

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');

% load LFP Information
[analogChannelsStored,~,timeVals] = loadlfpInfo(folderLFP);
[neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes);

% Get Combinations
[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);

% Get properties of the Stimulus
%stimResults = loadStimResults(folderExtract);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
% fonts
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.25; panelStartHeight = 0.7;
staticPanelWidth = 0.2; staticStartPos = 0.1;
dynamicPanelWidth = 0.2; dynamicStartPos = 0.3;
timingPanelWidth = 0.2; timingStartPos = 0.5;
plotOptionsPanelWidth = 0.2; plotOptionsStartPos = 0.7;
backgroundColor = 'w';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Static Panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%staticTitle = [subjectName '_' expDate '_' protocolName];
% if 0 % don't plot the static panel 
% hStaticPanel = uipanel('Title','Information','fontSize', fontSizeLarge, ...
%     'Unit','Normalized','Position',[staticStartPos panelStartHeight staticPanelWidth panelHeight]);
% 
% staticText = [{ '   '};
%     {['Monkey Name: ' subjectName]}; ...
%     {['Date: ' expDate]}; ...
%     {['Protocol Name: ' protocolName]}; ...
%     {'   '}
%     {['Orientation  (Deg): ' num2str(stimResults.orientation)]}; ...
%     {['Spatial Freq (CPD): ' num2str(stimResults.spatialFrequency)]}; ...
%     {['Eccentricity (Deg): ' num2str(stimResults.eccentricity)]}; ...
%     {['Polar angle  (Deg): ' num2str(stimResults.polarAngle)]}; ...
%     {['Sigma        (Deg): ' num2str(stimResults.sigma)]}; ...
%     {['Radius       (Deg): ' num2str(stimResults.radius)]}; ...
%     ];
% 
% tStaticText = uicontrol('Parent',hStaticPanel,'Unit','Normalized', ...
%     'Position',[0 0 1 1], 'Style','text','String',staticText,'FontSize',fontSizeSmall);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.06; dynamicGap=0.015; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Azimuth
azimuthString = getStringFromValues(aValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Azimuth (Deg)','FontSize',fontSizeSmall);
hAzimuth = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',azimuthString,'FontSize',fontSizeSmall);

% Elevation
elevationString = getStringFromValues(eValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-2*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Elevation (Deg)','FontSize',fontSizeSmall);
hElevation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',elevationString,'FontSize',fontSizeSmall);

% Sigma
sigmaString = getStringFromValues(sValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-3*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Sigma (Deg)','FontSize',fontSizeSmall);
hSigma = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-3*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',sigmaString,'FontSize',fontSizeSmall);

% Spatial Frequency
spatialFreqString = getStringFromValues(fValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-4*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeSmall);
hSpatialFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-4*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',spatialFreqString,'FontSize',fontSizeSmall);

% Orientation
orientationString = getStringFromValues(oValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Orientation (Deg)','FontSize',fontSizeSmall);
hOrientation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',orientationString,'FontSize',fontSizeSmall);

% Contrast
if ~isempty(cValsUnique)
    contrastString = getStringFromValues(cValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-6*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Contrast (%)','FontSize',fontSizeSmall);
    hContrast = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-6*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',contrastString,'FontSize',fontSizeSmall);
end

% Temporal Freq
if ~isempty(tValsUnique)
    temporalFreqString = getStringFromValues(tValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-7*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeSmall);
    hTemporalFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-7*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',temporalFreqString,'FontSize',fontSizeSmall);
end

% Analysis Type
analysisTypeString = 'ERP|Firing Rate|FFT|delta FFT';
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-8*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-8*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.1; timingTextWidth = 0.5; timingBoxWidth = 0.25;
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

signalRange = [-0.1 0.5];
fftRange = [0 250];
baseline = [-0.2 0];
stimPeriod = [0.2 0.4];

% Signal Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Parameter','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Min','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth+timingBoxWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Max','FontSize',fontSizeMedium);

% Stim Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim Range (s)','FontSize',fontSizeSmall);
hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(1)),'FontSize',fontSizeSmall);
hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(2)),'FontSize',fontSizeSmall);


% FFT Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-5*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','FFT Range (Hz)','FontSize',fontSizeSmall);
hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(1)),'FontSize',fontSizeSmall);
hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(2)),'FontSize',fontSizeSmall);

% Baseline
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Basline (s)','FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),'FontSize',fontSizeSmall);

% Stim Period
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-7*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim period (s)','FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

% Y Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-8*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Y Range','FontSize',fontSizeSmall);
hYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','0','FontSize',fontSizeSmall);
hYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','1','FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.1;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[plotOptionsStartPos panelStartHeight plotOptionsPanelWidth panelHeight]);

% Button for Plotting
[colorString, colorNames] = getColorString;
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Color','FontSize',fontSizeSmall);
hChooseColor = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 1-plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',colorString,'FontSize',fontSizeSmall);

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 4*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

hHoldOn = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium, ...
    'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale Y','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleY_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale X','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleData_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show electrode array and bad channels
electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];

% Get Bad channels
impedanceFileName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,'impedanceValues.mat');
if exist(impedanceFileName,'file')
    impedanceValues = getImpedanceValues(impedanceFileName);
    badChannels = find(impedanceValues>2500);
    highImpChannels = find(impedanceValues>1500);
else
    disp('Could not find impedance values');
    badChannels=[];
    highImpChannels=[];
end
     
hElectrodes = showElectrodeLocations(electrodeGridPos,highImpChannels,'m',[],0,0,gridType,subjectName,gridLayout);
              showElectrodeLocations(electrodeGridPos,badChannels,'r',[],1,0,gridType,subjectName,gridLayout);

% Get main plot and message handles
[~,~,electrodeArray] = electrodePositionOnGrid(1,gridType,subjectName,gridLayout);
[numRows,numCols] = size(electrodeArray);
plotHandles = getPlotHandles(numRows,numCols);

hMessage = uicontrol('Unit','Normalized','Position',[0 0.975 1 0.025],...
    'Style','text','String',[subjectName expDate protocolName],'FontSize',fontSizeLarge);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)
        a=get(hAzimuth,'val');
        e=get(hElevation,'val');
        s=get(hSigma,'val');
        f=get(hSpatialFreq,'val');
        o=get(hOrientation,'val');
        
        if ~isempty(tValsUnique)
            c=get(hContrast,'val');
            t=get(hTemporalFreq,'val');
            goodPos = parameterCombinations{a,e,s,f,o,c,t};
        elseif ~isempty(cValsUnique)
            c=get(hContrast,'val');
            goodPos = parameterCombinations{a,e,s,f,o,c};
        else
            goodPos = parameterCombinations{a,e,s,f,o};
        end
        
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));

        blRange = [str2double(get(hBaselineMin,'String')) str2double(get(hBaselineMax,'String'))];
        stRange = [str2double(get(hStimPeriodMin,'String')) str2double(get(hStimPeriodMax,'String'))];

        set(hMessage,'String',[num2str(length(goodPos)) ' stimuli found' ]);
        
        if ~isempty(goodPos)
            if analysisType == 2
                holdOnState = get(hHoldOn,'val');
                channelsStored = neuralChannelsStored;
                [baselineFiringRate,stimulusFiringRate] = plotSpikeData(plotHandles,channelsStored,goodPos,folderSpikes,...
                    timeVals,plotColor,SourceUnitID,holdOnState,blRange,stRange,gridType,subjectName);
                responsiveElectrodes = (channelsStored(stimulusFiringRate>=5));
                inhibitedElectrodes = (channelsStored(intersect(find(baselineFiringRate>=5),find(stimulusFiringRate<=5))));
                disp(['responsive: ' num2str(responsiveElectrodes)]);
                disp(['inhibited : ' num2str(inhibitedElectrodes)]);
                showElectrodeLocations(electrodeGridPos,responsiveElectrodes,'b',hElectrodes,1,0,gridType,subjectName,gridLayout);
                showElectrodeLocations(electrodeGridPos,inhibitedElectrodes,'g',hElectrodes,1,0,gridType,subjectName,gridLayout);
            
            else
                channelsStored = analogChannelsStored;
                plotLFPData(plotHandles,channelsStored,goodPos,folderLFP,...
                analysisType,timeVals,plotColor,blRange,stRange,gridType,subjectName,gridLayout);
            end
            
            if analysisType<=2 % ERP or spikes
                xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
            else
                xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
            end
            
            yRange = getYLims(plotHandles,channelsStored,gridType,subjectName,gridLayout);
            set(hYMin,'String',num2str(yRange(1))); set(hYMax,'String',num2str(yRange(2)));
            rescaleData(plotHandles,channelsStored,[xRange yRange],gridType,subjectName,gridLayout);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType == 2
            channelsStored = neuralChannelsStored;
        else
            channelsStored = analogChannelsStored;
        end
        
        if analysisType<=2 % ERP or spikes
            xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        else
            xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        end

        yRange = [str2double(get(hYMin,'String')) str2double(get(hYMax,'String'))];
        rescaleData(plotHandles,channelsStored,[xRange yRange],gridType,subjectName,gridLayout);

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType == 2
            channelsStored = neuralChannelsStored;
        else
            channelsStored = analogChannelsStored;
        end
        
        if analysisType<=2 % ERP or spikes
            xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        else
            xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        end

        yRange = getYLims(plotHandles,channelsStored,gridType,subjectName,gridLayout);
        rescaleData(plotHandles,channelsStored,[xRange yRange],gridType,subjectName,gridLayout);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');

        [numRow,numCol] = size(plotHandles);

        if holdOnState
            for i=1:numRow
                for j=1:numCol
                    set(plotHandles(i,j),'Nextplot','add');
                end
            end
        else
            for i=1:numRow
                for j=1:numCol
                    set(plotHandles(i,j),'Nextplot','replace');
                end
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        [numRow,numCol] = size(plotHandles);
        for i=1:numRow
            for j=1:numCol
                cla(plotHandles(i,j));
            end
        end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that plots the data
function plotLFPData(plotHandles, channelsStored, goodPos, ...
    folderData, analysisType, timeVals, plotColor,blRange,stRange,gridType,subjectName,gridLayout)

if isempty(goodPos)
    disp('No entries for this combination..')
else
    
    removeAvgRef = 0;
    
    if removeAvgRef
        disp('Removing average reference'); %#ok<*UNRCH>
        load([folderData 'avgRef']);
        avgRef = analogData;
    end
    
    if analysisType>2 % FFT
        Fs = round(1/(timeVals(2)-timeVals(1)));
        blPos = find(timeVals>=blRange(1),1)+ (1:diff(blRange)*Fs);
        stPos = find(timeVals>=stRange(1),1)+ (1:diff(stRange)*Fs);

        xsBL = 0:1/(diff(blRange)):Fs-1/(diff(blRange));
        xsST = 0:1/(diff(stRange)):Fs-1/(diff(stRange));
    end

    for i=1:length(channelsStored)
        
        channelNum = channelsStored(i);
        disp(['Plotting electrode ' num2str(channelNum)]);

        % get position
        [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName,gridLayout);

        if analysisType == 1        % compute ERP
            clear signal analogData
            load(fullfile(folderData ,['elec' num2str(channelNum)]));
            if removeAvgRef
                analogData = analogData-avgRef;
            end
            erp = mean(analogData(goodPos,:),1); %#ok<*NODEF>

            %Plot
            plot(plotHandles(row,column),timeVals,erp,'color',plotColor);

        elseif analysisType == 2    % compute Firing rates
            disp('Use plotSpikeData instead of plotLFPData...');
        else
            clear signal analogData
            load(fullfile(folderData,['elec' num2str(channelNum)]));
            if removeAvgRef
                analogData = analogData-avgRef;
            end
            fftBL = abs(fft(analogData(goodPos,blPos),[],2)); 
            fftST = abs(fft(analogData(goodPos,stPos),[],2));

            if analysisType == 3
                plot(plotHandles(row,column),xsBL,log10(mean(fftBL)),'g');
                set(plotHandles(row,column),'Nextplot','add');
                plot(plotHandles(row,column),xsST,log10(mean(fftST)),'k');
                set(plotHandles(row,column),'Nextplot','replace');
            end

            if analysisType == 4
                if xsBL == xsST %#ok<*BDSCI>
                    plot(plotHandles(row,column),xsBL,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor);
                    set(plotHandles(row,column),'Nextplot','add');
                    plot(plotHandles(row,column),xsBL,zeros(1,length(xsBL)),'color','k');
                    set(plotHandles(row,column),'Nextplot','replace');
                else
                    disp('Choose same baseline and stimulus periods..');
                end
            end
        end
    end
end
end

function [baselineFiringRate,stimulusFiringRate] = plotSpikeData(plotHandles,channelsStored,goodPos, ...
    folderData, timeVals, plotColor, SourceUnitID,holdOnState,blRange,stRange,gridType,subjectName)

unitColors = ['r','m','y','c','g'];
binWidthMS = 10;

if isempty(goodPos)
    disp('No entries for this combination..')
else
    numChannels = length(channelsStored);
    baselineFiringRate = zeros(1,numChannels);
    stimulusFiringRate = zeros(1,numChannels);
    
    for i=1:length(channelsStored)
        channelNum = channelsStored(i);
        disp(channelNum)

        % get position
        [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName);

        clear neuralInfo spikeData
        load(fullfile(folderData,['elec' num2str(channelNum) '_SID' num2str(SourceUnitID(i))]));
        [psthVals,xs] = getPSTH(spikeData(goodPos),binWidthMS,[timeVals(1) timeVals(end)]);
        
        % Compute the mean firing rates
        blPos = find(xs>=blRange(1),1)+ (1:(diff(blRange))/(binWidthMS/1000));
        stPos = find(xs>=stRange(1),1)+ (1:(diff(stRange))/(binWidthMS/1000));
        
        baselineFiringRate(i) = mean(psthVals(blPos));
        stimulusFiringRate(i) = mean(psthVals(stPos));
        
        if SourceUnitID(i)==0
            plot(plotHandles(row,column),xs,psthVals,'color',plotColor);
        elseif SourceUnitID(i)> 5
            disp('Only plotting upto 6 single units per electrode...')
        else
            plot(plotHandles(row,column),xs,smooth(psthVals),'color',unitColors(SourceUnitID(i)));
        end
        
        if (i<length(channelsStored))
            if channelsStored(i) == channelsStored(i+1)
                disp('hold on...')
                set(plotHandles(row,column),'Nextplot','add');
            else
                if ~holdOnState
                    set(plotHandles(row,column),'Nextplot','replace');
                end
            end
        end
        
    end
end
end   

function yRange = getYLims(plotHandles,channelsStored,gridType,subjectName,gridLayout)

% Initialize
yMin = inf;
yMax = -inf;

for i=1:length(channelsStored)
    channelNum = channelsStored(i);
    % get position
    [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName,gridLayout);
    
    axis(plotHandles(row,column),'tight');
    tmpAxisVals = axis(plotHandles(row,column));
    if tmpAxisVals(3) < yMin
        yMin = tmpAxisVals(3);
    end
    if tmpAxisVals(4) > yMax
        yMax = tmpAxisVals(4);
    end
end
yRange = [yMin yMax];
end
function rescaleData(plotHandles,channelsStored,axisLims,gridType,subjectName,gridLayout)

[numRows,numCols] = size(plotHandles);
labelSize=12;
for i=1:length(channelsStored)
    channelNum = channelsStored(i);
    
    % get position
   [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName,gridLayout);
    
    axis(plotHandles(row,column),axisLims);
    if (row==numRows && rem(column,2)==rem(numCols+1,2))
        if column~=1
            set(plotHandles(row,column),'YTickLabel',[],'fontSize',labelSize);
        end
    elseif (rem(row,2)==0 && column==1)
        set(plotHandles(row,column),'XTickLabel',[],'fontSize',labelSize);
    else
        set(plotHandles(row,column),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
    end
end

% Remove Labels on the four corners
set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%
%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Data
function [analogChannelsStored,goodStimPos,timeVals,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end
function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
fileName = fullfile(folderSpikes,'spikeInfo.mat');
if exist(fileName,'file')
    load(fileName);
else
    neuralChannelsStored=[];
    SourceUnitID=[];
end
end
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('sValsUnique','var');    sValsUnique = rValsUnique/3;         end
if ~exist('cValsUnique','var');    cValsUnique=[];                      end
if ~exist('tValsUnique','var');    tValsUnique=[];                      end
end
% function stimResults = loadStimResults(folderExtract)
% load ([folderExtract 'stimResults']);
% end
function impedanceValues = getImpedanceValues(fileName)
load(fileName);
end