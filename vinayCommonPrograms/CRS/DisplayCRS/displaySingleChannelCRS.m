% Displays data from a single electrode, bipolar signals
% Vinay - modified from displaySingleChannelGRFv3.m for CRS protocol
% 09 August 2014

function displaySingleChannelCRS(subjectName,expDate,protocolName,folderSourceString,gridType,loadProtocolNumber)

if ~exist('folderSourceString','var')  folderSourceString='/media/store/';        end
if ~exist('gridType','var')            gridType='EEG';                  end
if ~exist('loadProtocolNumber','var')  loadProtocolNumber = 11;         end % Vinay - anything greater than 10 will read the protocolNumber from the saved data


% Get folders
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');
%-----

% load LFP Information
[analogChannelsStored,timeVals,~,analogInputNums] = loadlfpInfo(folderLFP);
[neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes);

% Get Combinations
[~,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);

% Get properties of the Stimulus
% stimResults = loadStimResults(folderExtract);

% Vinay - read the variable loadProtocolNumber as protocolNumber if it is
% less than 11. Otherwise it will be loaded from stimResults. This is so
% that one can directly pass the protocolNumber in case it is not recorded
% in stimResults
if (loadProtocolNumber < 11)
    protocolNumber = loadProtocolNumber;
else
    protocolNumber = getProtocolNumber(folderExtract);
end

% Vinay - get the particular gabors that were displayed in the protocol.
gaborsDisplayed = getGaborsDisplayed(folderExtract);
disp(['Gabors displayed in this protocol:' num2str(gaborsDisplayed)]);

% Vinay - if all the gabors were not displayed then assign some default
% value to the null gabor's parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
% fonts
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; fontSizeTiny = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.34; panelStartHeight = 0.61;
staticPanelWidth = 0.22; staticStartPos = 0.025; % [Vinay] - changed width 0.25 to 0.22
dynamicPanelWidth = 0.23; dynamicStartPos = 0.245; % [Vinay] - changed width 0.25 to 0.23
timingPanelWidth = 0.20; timingStartPos = 0.475; % [Vinay] - changed width 0.25 to 0.20
plotOptionsPanelWidth = 0.15; plotOptionsStartPos = 0.675; % [Vinay] - changed width 0.2 to 0.15
selectParamPanelWidth = 0.15; selectParamStartPos = 0.825; % [Vinay] - adding a panel to select parameters
backgroundColor = 'w';

% Vinay - define a flag for notching the line noise
notchData = 0;
% [Vinay] - define a flag to take bipolar (i.e. diff) signals: V1-V2
useBipolar = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.06; dynamicGap=0.015; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Analog channel
[analogChannelStringList,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Analog Channel (1:2)','FontSize',fontSizeTiny);
hAnalogChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) (0.5*(1-dynamicTextWidth)) dynamicHeight], ...
    'Style','popup','String',analogChannelStringList,'FontSize',fontSizeTiny);
% [Vinay] - adding analog channel 2 for chosing bipolar signals V1-V2 
hAnalogChannel2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(dynamicTextWidth+(0.5*(1-dynamicTextWidth))) 1-(dynamicHeight+dynamicGap) (0.5*(1-dynamicTextWidth)) dynamicHeight], ...
    'Style','popup','String',(['none|' analogChannelStringList]),'FontSize',fontSizeTiny,'Callback',{@bipolar_Callback});

% Neural channel
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-2*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Neural Channel','FontSize',fontSizeTiny);
    
if ~isempty(neuralChannelsStored)
    neuralChannelString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs);
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',neuralChannelString,'FontSize',fontSizeTiny);
else
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position', [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Not found','FontSize',fontSizeTiny);
end

%[Vinay] - create 3 columns corresponding to the parameters of the 3 gabors
%-Centre, Ring, Surround. Textwidth = 0.4, C,R,S popout - each 0.2 width

labelWidth = 0.4; crsWidth = 0.2; % [Vinay] - labelWidth + 3(crsWidth) = 1

uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[labelWidth 1-3*(dynamicHeight+dynamicGap) (3*crsWidth) dynamicHeight],...
    'Style','text','String','Surround        Ring        Centre','FontSize',fontSizeTiny);

for gaborNum = 1:3
    
    % Azimuth
    azimuthString = getStringFromValues(aValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-4*(dynamicHeight+dynamicGap) labelWidth dynamicHeight],...
    'Style','text','String','Azimuth (Deg)','FontSize',fontSizeTiny);
    hAzimuth(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-4*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',azimuthString,'FontSize',fontSizeTiny);

    % Elevation
    elevationString = getStringFromValues(eValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-5*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Elevation (Deg)','FontSize',fontSizeTiny);
    hElevation(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-5*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',elevationString,'FontSize',fontSizeTiny);
    
    % Sigma
    sigmaString = getStringFromValues(sValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-6*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Sigma (Deg)','FontSize',fontSizeTiny);
    hSigma(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-6*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',sigmaString,'FontSize',fontSizeTiny);

    % Spatial Frequency
    spatialFreqString = getStringFromValues(fValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-7*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeTiny);
    hSpatialFreq(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-7*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',spatialFreqString,'FontSize',fontSizeTiny);

    % Orientation
    orientationString = getStringFromValues(oValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-8*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Orientation (Deg)','FontSize',fontSizeTiny);
    hOrientation(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-8*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',orientationString,'FontSize',fontSizeTiny);

    % Contrast
    contrastString = getStringFromValues(cValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-9*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Contrast (%)','FontSize',fontSizeTiny);
    hContrast(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-9*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',contrastString,'FontSize',fontSizeTiny);

    % Temporal Frequency
    temporalFreqString = getStringFromValues(tValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-10*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeTiny);
    hTemporalFreq(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-10*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',temporalFreqString,'FontSize',fontSizeTiny);
    
    % [Vinay] - adding radius and spatial phase cells
    % radius
    radiusString = getStringFromValues(rValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-11*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Radius (deg)','FontSize',fontSizeTiny);
    hRadius(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-11*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',radiusString,'FontSize',fontSizeTiny);

    % spatial phase
    spatialPhaseString = getStringFromValues(pValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-12*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Spatial Phase (Deg)','FontSize',fontSizeTiny);
    hSpatialPhase(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-12*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',spatialPhaseString,'FontSize',fontSizeTiny);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.1; timingTextWidth = 0.5; timingBoxWidth = 0.20; % [Vinay] - changed width from 0.25 to 0.20 
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

% Default values
signalRange = [-0.8 0.9];
fftRange = [0 250];
baseline = [-0.5 0];
stimPeriod = [0.2 0.7];

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

% STA length
staLen = [-0.05 0.05]; 
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-9*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','STA Len(i)(s)','FontSize',fontSizeSmall);
hSTAMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(staLen(1)),'FontSize',fontSizeSmall);
hSTAMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(staLen(2)),'FontSize',fontSizeSmall);
hRemoveMeanSTA = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0 1-10*timingHeight 1 timingHeight], ...
    'Style','togglebutton','String','remove mean STA','FontSize',fontSizeMedium);

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

% Analysis Type
analysisTypeString = 'ERP|Firing Rate|Raster|FFT|delta FFT|STA';
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-3*plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.6 1-3*plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall);


% Vinay - adding a toggle button to use notched data
hNotchData = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 5*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','raw/lineNoiseRemoved','FontSize',fontSizeMedium, ...
    'Callback',{@notch_Callback});

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
% [Vinay] - adding a panel to select the parameters for the plots

selectParamHeight = 0.1; selectParamWidth = 0.5; selectParamBoxWidth = 0.20; % [Vinay] - changed width from 0.25 to 0.20
textWidth = 0.4; selGWidth = 0.2; selPWidth = 0.4;
hSelectParamPanel = uipanel('Title','Select Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[selectParamStartPos panelStartHeight selectParamPanelWidth panelHeight]);


parametersString = 'NA|azimuth|elevation|sigma|spatialFreq|orientation|contrast|temporalFreq|radius|spatialPhase';
gaborString = 'NA|C|R|S';

%parameter1
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-selectParamHeight textWidth selectParamHeight], ...
    'Style','text','String','param1','FontSize',fontSizeSmall);
hGabor1 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [textWidth 1-selectParamHeight selGWidth selectParamHeight], ...
    'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});
hParam1 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+selGWidth) 1-selectParamHeight selPWidth selectParamHeight], ...
    'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});


%parameter2
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-2*selectParamHeight textWidth selectParamHeight], ...
    'Style','text','String','param2','FontSize',fontSizeSmall);
hGabor2 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [textWidth 1-2*selectParamHeight selGWidth selectParamHeight], ...
    'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});
hParam2 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+selGWidth) 1-2*selectParamHeight selPWidth selectParamHeight], ...
    'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});


%parameter3
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-3*selectParamHeight textWidth selectParamHeight], ...
    'Style','text','String','param3','FontSize',fontSizeSmall);
hGabor3 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [textWidth 1-3*selectParamHeight selGWidth selectParamHeight], ...
    'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});
hParam3 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+selGWidth) 1-3*selectParamHeight selPWidth selectParamHeight], ...
    'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});


%parameter4
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-4*selectParamHeight textWidth selectParamHeight], ...
    'Style','text','String','param4','FontSize',fontSizeSmall);
hGabor4 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [textWidth 1-4*selectParamHeight selGWidth selectParamHeight], ...
    'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});
hParam4 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+selGWidth) 1-4*selectParamHeight selPWidth selectParamHeight], ...
    'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});


%parameter5
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-5*selectParamHeight textWidth selectParamHeight], ...
    'Style','text','String','param5','FontSize',fontSizeSmall);
hGabor5 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [textWidth 1-5*selectParamHeight selGWidth selectParamHeight], ...
    'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});
hParam5 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+selGWidth) 1-5*selectParamHeight selPWidth selectParamHeight], ...
    'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});


%parameter6
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-6*selectParamHeight textWidth selectParamHeight], ...
    'Style','text','String','param6','FontSize',fontSizeSmall);
hGabor6 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [textWidth 1-6*selectParamHeight selGWidth selectParamHeight], ...
    'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});
hParam6 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+selGWidth) 1-6*selectParamHeight selPWidth selectParamHeight], ...
    'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});


%Conjunction plot parameter and value. This parameter will be fixed to a
%value and the rest of the plots are taken by averaging across only this
%value of this parameter. Eg. All plots at contrast 50% of centre gabor
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-7*selectParamHeight textWidth selectParamHeight], ...
    'Style','text','String','pivot','FontSize',fontSizeSmall);
hGabor7 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [textWidth 1-7*selectParamHeight selGWidth selectParamHeight], ...
    'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});
hParam7 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+selGWidth) 1-7*selectParamHeight selPWidth selectParamHeight], ...
    'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetParams_Callback});


% check box - 'use default'. This will load default parameters as per the
% specific protocol being used
useDefParams = 0; % by default don't use the default parameters
hUseDefParams = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.1 1-9*selectParamHeight 0.8 selectParamHeight], ...
    'Style','checkbox','String','Use default parameters','FontSize',fontSizeTiny, 'Callback',{@useDefParams_Callback});

% display the specific protocol being used
protocolNameString = {'NoneProtocol';'RingProtocol';'ContrastRingProtocol';'DualContrastProtocol'; ...
    'DualOrientationProtocol';'DualPhaseProtocol';'OrientationRingProtocol';'PhaseRingProtocol'; ...
    'Drifting Ring Protocol';'CrossOrientationProtocol';'AnnulusRingProtocol'};
% protocolNumber = getProtocolNumber(folderExtract);

uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0.1 1-10*selectParamHeight 0.8 selectParamHeight], ...
    'Style','text','String',protocolNameString(protocolNumber+1),'FontSize',fontSizeTiny,'FontWeight','bold');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get plots and message handles

% Get electrode array information
electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];
hElectrodes = showElectrodeLocations(electrodeGridPos,analogChannelsStored(get(hAnalogChannel,'val')), ...
    colorNames(get(hChooseColor,'val')),[],1,0,gridType);

% Make plot for RFMap, centerRFMap and main Map
if length(aValsUnique)>=5
    mapRatio = 2/3; % this sets the relative ratio of the mapping plots versus orientation plots
else
    mapRatio = 1/2;
end

startXPos = staticStartPos; endXPos = 0.95; startYPos = 0.05; mainRFHeight = 0.55; centerGap = 0.05;
mainRFWidth = mapRatio*(endXPos-startXPos-centerGap);
otherPlotsWidth = (1-mapRatio)*(endXPos-startXPos-centerGap);

% RF and centerRF
RFMapPos = [endXPos-mainRFWidth startYPos+(3/4)*mainRFHeight mainRFWidth/2 (1/4)*mainRFHeight];
hRFMapPlot = subplot('Position',RFMapPos,'XTickLabel',[],'YTickLabel',[],'box','on');
centerRFMapPos = [endXPos-mainRFWidth+mainRFWidth/2 startYPos+(3/4)*mainRFHeight mainRFWidth/2 (1/4)*mainRFHeight];
hcenterRFMapPlot = subplot('Position',centerRFMapPos,'XTickLabel',[],'YTickLabel',[],'box','on');

% Main plot handles 
numRows = length(eValsUnique{1}); numCols = length(aValsUnique{1}); % [Vinay] -since C,R,S are concentric, 
% the azi and ele of surround are sufficient to cover the cases, hence e/aValsUnique{1}
gridPos=[endXPos-mainRFWidth startYPos mainRFWidth 0.65*mainRFHeight]; gap = 0.002;
plotHandles = getPlotHandles(numRows,numCols,gridPos,gap);

uicontrol('Unit','Normalized','Position',[0 0.975 1 0.025],...
    'Style','text','String',[subjectName expDate protocolName],'FontSize',fontSizeLarge);

% Other functions

% Remaining Grid size
remainingWidth = otherPlotsWidth;
remainingHeight= mainRFHeight;

otherGapSize = 0.04;
        
otherHeight = (remainingHeight-5*otherGapSize)/6;

param1Grid = [startXPos startYPos                               remainingWidth otherHeight];
param2Grid = [startXPos startYPos+ (otherHeight+otherGapSize)   remainingWidth otherHeight];
param3Grid = [startXPos startYPos+ 2*(otherHeight+otherGapSize) remainingWidth otherHeight];
param4Grid = [startXPos startYPos+ 3*(otherHeight+otherGapSize) remainingWidth otherHeight];
param5Grid = [startXPos startYPos+ 4*(otherHeight+otherGapSize) remainingWidth otherHeight];
param6Grid = [startXPos startYPos+ 5*(otherHeight+otherGapSize) remainingWidth otherHeight];

% Plot handles
% [Vinay] - get the plot handles based on the parameters selected. This
% should change everytime one selects any gabor and parameter. Therefore
% this is now done in the callbacks. Default is done below
% hParam1Plot = getPlotHandles(1,length(tValsUnique),param1Grid,0.002);
% hParam2Plot = getPlotHandles(1,length(cValsUnique),param2Grid,0.002);
% hParam3Plot = getPlotHandles(1,length(oValsUnique),param4Grid,0.002);
% hParam4Plot = getPlotHandles(1,length(fValsUnique),param3Grid,0.002);
% hParam5Plot = getPlotHandles(1,length(sValsUnique),param5Grid,0.002);
% hParam5Plot = getPlotHandles(1,length(sValsUnique),param6Grid,0.002);

% [Vinay] - load default parameters
set(hGabor1,'val',4); % 4 => S
set(hGabor2,'val',4);
set(hGabor3,'val',3); % 3 => R
set(hGabor4,'val',3);
set(hGabor5,'val',2); % 2 => C
set(hGabor6,'val',2);
set(hGabor7,'val',1); % 1 => NA

set(hParam1,'val',7); % 7 => cont (surr)
set(hParam2,'val',6); % 6 => ori
set(hParam3,'val',7); % 7 => cont (ring)
set(hParam4,'val',9); % 9 => rad 
set(hParam5,'val',7); % 7 => cont (centre)
set(hParam6,'val',6); % 6 => ori
set(hParam7,'val',1); % 1 => NA

adjFactor = 5;
gabor1 = adjFactor - get(hGabor1,'val'); % If sel is C, then value = 2 (NA = 1, C = 2, R = 3, S = 4) => gabor1 = 5 - 2 = 3. If R then 2, if S then 1
gabor2 = adjFactor - get(hGabor2,'val');
gabor3 = adjFactor - get(hGabor3,'val');
gabor4 = adjFactor - get(hGabor4,'val');
gabor5 = adjFactor - get(hGabor5,'val');
gabor6 = adjFactor - get(hGabor6,'val');
gabor7 = adjFactor - get(hGabor7,'val');

param1 = get(hParam1, 'val');
param2 = get(hParam2, 'val');
param3 = get(hParam3, 'val');
param4 = get(hParam4, 'val');
param5 = get(hParam5, 'val');
param6 = get(hParam6, 'val');
param7 = get(hParam7, 'val');

hParam1Plot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1Grid,0.002);
hParam2Plot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2Grid,0.002);
hParam3Plot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3Grid,0.002);
hParam4Plot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4Grid,0.002);
hParam5Plot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5Grid,0.002);
hParam6Plot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6Grid,0.002);


%**************************************************************************
% [Vinay] - TIME-FREQUENCY ANALYSIS
% Plot the time-frequency distributions on a separate figure
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hTFfig = figure(2);

% Define variables/flags required in TF analysis with their default values
plotStyle = 3; % pcolor, imagesc, raw
spectrumType = 1; % raw, difference
cmin = -1; cmax = 3; % caxis limits
tfMethod = 1; % MTM, MP

% Default MTM params
mtmParams.Fs = 2000;
mtmParams.tapers=[1 1]; % [1 1] case is simply stft with dpss window
mtmParams.trialave=1;
mtmParams.err=0;
mtmParams.pad=0;

movingWin = [0.5 0.01];

% Default MP parameters
numAtomsMP = 100;

%__________________________________________________________________________
% The tf plots panel for the 6 chosen parameters
% ---------------tfPlotsPanel----------------------------------------------
%--------------------------------------------------------------------------

xGap = 0.01; yGap = 0.01;
tfPlotsHeight = 0.14; tfPlotsGap =0.02;
tfPanelWidth = 0.6; tfPanelHeight = 1;
hTFPlotsPanel = uipanel('Title','TF Plots','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[xGap yGap tfPanelWidth-2*xGap tfPanelHeight-2*yGap]);

startXPos = 2*xGap; startYPos = 2*yGap; % start from left bottom corner
tfPanelRowWidth = tfPanelWidth-4*xGap;

param1TFGrid = [startXPos startYPos                               tfPanelRowWidth tfPlotsHeight];
param2TFGrid = [startXPos startYPos+ (tfPlotsHeight+tfPlotsGap)   tfPanelRowWidth tfPlotsHeight];
param3TFGrid = [startXPos startYPos+ 2*(tfPlotsHeight+tfPlotsGap) tfPanelRowWidth tfPlotsHeight];
param4TFGrid = [startXPos startYPos+ 3*(tfPlotsHeight+tfPlotsGap) tfPanelRowWidth tfPlotsHeight];
param5TFGrid = [startXPos startYPos+ 4*(tfPlotsHeight+tfPlotsGap) tfPanelRowWidth tfPlotsHeight];
param6TFGrid = [startXPos startYPos+ 5*(tfPlotsHeight+tfPlotsGap) tfPanelRowWidth tfPlotsHeight];


hParam1TFPlot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1TFGrid,0.002);
hParam2TFPlot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2TFGrid,0.002);
hParam3TFPlot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3TFGrid,0.002);
hParam4TFPlot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4TFGrid,0.002);
hParam5TFPlot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5TFGrid,0.002);
hParam6TFPlot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6TFGrid,0.002);


%__________________________________________________________________________
% Parameters selection panels for TF methods
% ---------------tfParamsPanel----------------------------------------------
%--------------------------------------------------------------------------

xGap = 0.002;
tfParamsHeight = tfPlotsHeight; tfParamsGap = tfPlotsGap;
tfParamsWidth = 1 - (tfPanelWidth+xGap); tfParamsPanelHeight = 2*(tfParamsHeight+tfParamsGap);
tfParamsPanelxStart = (tfPanelWidth+xGap);
tfParamsPanelyStart = (1-tfParamsPanelHeight);

% hTFParamsPanel = uipanel('Title','TF Parameters','fontSize', fontSizeLarge, ...
%     'Unit','Normalized','Position',[tfParamsPanelxStart tfParamsPanelyStart tfParamsWidth tfParamsPanelHeight]);

textWidth = 0.24;
textHeight = 0.18;

% ====================TF Params Panel=====================================
hTFParamsPanel = uipanel('Title','MTM, MP Parameters','fontSize', fontSizeMedium, ...
    'Unit','Normalized','Position',[(tfParamsPanelxStart) (tfParamsPanelyStart) tfParamsWidth/2 tfParamsPanelHeight]);

% Parameters to divide the panel into rows and columns
xPanel1 = 0; xPanel2 = 0.5;

% Fs, sampling frequency
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-textHeight) textWidth textHeight], ...
    'Style','text','String','Fs','FontSize',fontSizeSmall);

hMTMFs = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-textHeight) textWidth textHeight], ...
    'Style','edit','String',mtmParams.Fs,'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize with default Fs = 2000

%------------------MTM parameters-----------------------
% Tapers TW
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-2*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','text','String','TW','FontSize',fontSizeSmall);

hMTMTapersTW = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-2*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',mtmParams.tapers(1),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize TW = 1

% Tapers 'k'
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','k','FontSize',fontSizeSmall);

hMTMTapersK = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',mtmParams.tapers(2),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize k = 1

% Window length
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-4*textHeight-3*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','wLen','FontSize',fontSizeSmall);

hMTMwLen = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-4*textHeight-3*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',movingWin(1),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize wLen = 0.1 s

% Window translation step
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-5*textHeight-4*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','wStep','FontSize',fontSizeSmall);
hMTMwStep = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-5*textHeight-4*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',movingWin(2),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize wStep = 0.01 s

%----MP parameters-----
% number of atoms for reconstruction
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-textHeight) textWidth textHeight], ...
    'Style','text','String','numAtoms MP','FontSize',fontSizeSmall);
hMPnumAtoms = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel2+textWidth) (1-textHeight) textWidth textHeight], ...
    'Style','edit','String',numAtomsMP,'FontSize',fontSizeTiny,'Callback',{@resetMPParams_Callback}); % initialize numAtoms = 100


%========= Plotting settings ==============================
hTFPlotSettingsPanel = uipanel('Title','TF Plot settings','fontSize', fontSizeMedium, ...
    'Unit','Normalized','Position',[(tfParamsPanelxStart+(tfParamsWidth/2)+xGap) (tfParamsPanelyStart) ((tfParamsWidth/2)-xGap) tfParamsPanelHeight]);

% pcolor or imagesc or raw
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-textHeight) 2*textWidth textHeight], ...
    'Style','text','String','Plot style','FontSize',fontSizeSmall);
hTFPlotStyle = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-textHeight) 2*textWidth textHeight], ...
    'Style','pop','String','pcolor|imagesc|line','FontSize',fontSizeTiny,'Callback',{@resetTFSettings_Callback}); % initialize style = pcolor

% spectrum type: raw or difference (change from baseline)
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-2*textHeight-tfParamsGap) 2*textWidth textHeight], ...
    'Style','text','String','Spectrum type','FontSize',fontSizeSmall);
hTFSpectrumType = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-2*textHeight-tfParamsGap) 2*textWidth textHeight], ...
    'Style','pop','String','raw|difference','FontSize',fontSizeTiny,'Callback',{@resetTFSettings_Callback}); % initialize style = pcolor


% color axis settings
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','cmin','FontSize',fontSizeSmall);
hTFcmin = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',cmin,'FontSize',fontSizeTiny,'Callback',{@rescaleCaxis_Callback}); % initialize cmin = -1

uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','cmax','FontSize',fontSizeSmall);
hTFcmax = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel2+textWidth) (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',cmax,'FontSize',fontSizeTiny,'Callback',{@rescaleCaxis_Callback}); % initialize cmax = 3


% Method used
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-5*textHeight-4*tfParamsGap) 2*textWidth 1.5*textHeight], ...
    'Style','text','String','Method','FontSize',fontSizeMedium);
hTFPlotMethod = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-5*textHeight-4*tfParamsGap) 2*textWidth 1.5*textHeight], ...
    'Style','pop','String','MTM|MP','FontSize',fontSizeTiny,'Callback',{@resetTFSettings_Callback}); % initialize method = MTM


% --------- TF Plot callback button ---------------------------
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 textHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotTFData_Callback});

%==========================================================================
%-----------Main TF plot box---------------------
xGap = 0.005; yGap = 0.05;
numRows = length(eValsUnique{1}); numCols = length(aValsUnique{1}); % [Vinay] -since C,R,S are concentric, 
% the azi and ele of surround are sufficient to cover the cases, hence e/aValsUnique{1}
tfgridPos=[(tfPanelWidth+xGap) yGap (1-tfPanelWidth-4*xGap) (1-tfParamsPanelHeight-xGap-2*yGap)]; gap = 0.002;
tfplotHandles = getPlotHandles(numRows,numCols,tfgridPos,gap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)
        a=cell2mat(get(hAzimuth,'val'));
        e=cell2mat(get(hElevation,'val'));
        s=cell2mat(get(hSigma,'val'));
        f=cell2mat(get(hSpatialFreq,'val'));
        o=cell2mat(get(hOrientation,'val'));
        c=cell2mat(get(hContrast,'val'));
        t=cell2mat(get(hTemporalFreq,'val'));
        r=cell2mat(get(hRadius,'val'));
        p=cell2mat(get(hSpatialPhase,'val'));
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        BLMin = str2double(get(hBaselineMin,'String'));
        BLMax = str2double(get(hBaselineMax,'String'));
        STMin = str2double(get(hStimPeriodMin,'String'));
        STMax = str2double(get(hStimPeriodMax,'String'));
        STAMin = str2double(get(hSTAMin,'String'));
        STAMax = str2double(get(hSTAMax,'String'));
        holdOnState = get(hHoldOn,'val');
        removeMeanSTA = get(hRemoveMeanSTA,'val');
        notchData = get(hNotchData, 'val');
        
        
        % [Vinay] - get the plot handles here again
        hParam1Plot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1Grid,0.002);
        hParam2Plot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2Grid,0.002);
        hParam3Plot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3Grid,0.002);
        hParam4Plot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4Grid,0.002);
        hParam5Plot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5Grid,0.002);
        hParam6Plot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6Grid,0.002);


        if analysisType==6 % Spike triggered average
            analogChannelPos = get(hAnalogChannel,'val');
            analogChannelString = analogChannelStringArray{analogChannelPos};
            spikeChannelPos = get(hNeuralChannel,'val');
            spikeChannelNumber = neuralChannelsStored(spikeChannelPos);
            unitID = SourceUnitIDs(spikeChannelPos);
            
            plotColors{1} = 'g';
            plotColors{2} = 'k';
            plotSTA1Channel(plotHandles,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
                s,f,o,c,t,timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName,[STAMin STAMax],removeMeanSTA);
            
            % Write code for this
            %plotSTA1Parameter1Channel(hParam3Plot,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
            %    a,e,s,f,[],timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName);
            %plotSTA1Parameter1Channel(hParam4Plot,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
            %    a,e,s,[],o,timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName);
            
            if analogChannelPos<=length(analogChannelsStored)
                analogChannelNumber = analogChannelsStored(analogChannelPos);
            else
                analogChannelNumber = 0;
            end
            channelNumber = [analogChannelNumber spikeChannelNumber];
            
        elseif analysisType == 2 || analysisType == 3
            channelPos = get(hNeuralChannel,'val');
            channelNumber = neuralChannelsStored(channelPos);
            unitID = SourceUnitIDs(channelPos);
            plotSpikeData1Channel(plotHandles,channelNumber,s,f,o,c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
            plotSpikeData1Parameter1Channel(hParam1Plot,channelNumber,a,e,s,f,o,c,[],folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
            plotSpikeData1Parameter1Channel(hParam2Plot,channelNumber,a,e,s,f,o,[],t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
            plotSpikeData1Parameter1Channel(hParam3Plot,channelNumber,a,e,s,f,[],c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
            plotSpikeData1Parameter1Channel(hParam4Plot,channelNumber,a,e,s,[],o,c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
            plotSpikeData1Parameter1Channel(hParam5Plot,channelNumber,a,e,[],f,o,c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
        else
            analogChannelPos = get(hAnalogChannel,'val');
            analogChannelString = analogChannelStringArray{analogChannelPos};
            
            % [Vinay] - read the analog channel 2 string for bipolar case
            analogChannelPos2 = get(hAnalogChannel2,'val');
            analogChannelString2 = ('none');
            if analogChannelPos2 ~= 1
                analogChannelString2 = analogChannelStringArray{analogChannelPos2-1};
            end
            
            rfMapVals = plotLFPData1Channel(plotHandles,analogChannelString,analogChannelString2,s,f,o,c,t,r,p,folderLFP,...
                analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData,useBipolar); % Vinay - added notchData
            plotLFPData1Parameter1Channel(hParam1Plot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor1,param1,folderLFP,...
                analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar); % Vinay - added notchData
            plotLFPData1Parameter1Channel(hParam2Plot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor2,param2,folderLFP,...
                analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar); % Vinay - added notchData
            plotLFPData1Parameter1Channel(hParam3Plot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor3,param3,folderLFP,...
                analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar); % Vinay - added notchData
            plotLFPData1Parameter1Channel(hParam4Plot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor4,param4,folderLFP,...
                analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar); % Vinay - added notchData
            plotLFPData1Parameter1Channel(hParam5Plot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor5,param5,folderLFP,...
                analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar); % Vinay - added notchData
            plotLFPData1Parameter1Channel(hParam6Plot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor6,param6,folderLFP,...
                analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar); % Vinay - added notchData

            if analogChannelPos<=length(analogChannelsStored)
                channelNumber = analogChannelsStored(analogChannelPos);
            else
                channelNumber = 0;
            end
            
            if ~isempty(rfMapVals)
                if (length(aValsUnique{1})==1) || (length(eValsUnique{1})==1) % [Vinay] - using the azi and ele from S only since all are concentric
                    disp('Not enough data to plot RF center...')
                else
                    plotRFMaps(hRFMapPlot,hcenterRFMapPlot,rfMapVals,aValsUnique{1},eValsUnique{1},plotColor,holdOnState); % [Vinay] - using the azi and ele from S only since all are concentric
                end
            end
        end

        if analysisType<=3  % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        elseif analysisType <=5
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        else
            xMin = str2double(get(hSTAMin,'String'));
            xMax = str2double(get(hSTAMax,'String'));
        end

        rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
        rescaleData(hParam1Plot,xMin,xMax,getYLims(hParam1Plot));
        rescaleData(hParam2Plot,xMin,xMax,getYLims(hParam2Plot));
        rescaleData(hParam3Plot,xMin,xMax,getYLims(hParam3Plot));
        rescaleData(hParam4Plot,xMin,xMax,getYLims(hParam4Plot));
        rescaleData(hParam5Plot,xMin,xMax,getYLims(hParam5Plot));
        rescaleData(hParam6Plot,xMin,xMax,getYLims(hParam6Plot));
        showElectrodeLocations(electrodeGridPos,channelNumber,plotColor,hElectrodes,holdOnState,0,gridType);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');
        
        if analysisType<=3 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        else
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end

        yLims = [str2double(get(hYMin,'String')) str2double(get(hYMax,'String'))];
        rescaleData(plotHandles,xMin,xMax,yLims);
        rescaleData(hParam1Plot,xMin,xMax,yLims);
        rescaleData(hParam2Plot,xMin,xMax,yLims);
        rescaleData(hParam3Plot,xMin,xMax,yLims);
        rescaleData(hParam4Plot,xMin,xMax,yLims);
        rescaleData(hParam5Plot,xMin,xMax,yLims);
        rescaleData(hParam6Plot,xMin,xMax,yLims);
        
        % [Vinay] - rescale the TF plots
        figure(2);
        
        if plotStyle ~= 1 % pcolor, imagesc
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
            yLims = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
            rescaleData(hParam1TFPlot,xMin,xMax,yLims);
            rescaleData(hParam2TFPlot,xMin,xMax,yLims);
            rescaleData(hParam3TFPlot,xMin,xMax,yLims);
            rescaleData(hParam4TFPlot,xMin,xMax,yLims);
            rescaleData(hParam5TFPlot,xMin,xMax,yLims);
            rescaleData(hParam6TFPlot,xMin,xMax,yLims);

            for i = 1:size(tfplotHandles,1)
                for j = 1:size(tfplotHandles,2)
                    rescaleData(tfplotHandles(i,j),xMin,xMax,yLims);
                end
            end
            
        else % line
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
            rescaleData(hParam1TFPlot,xMin,xMax,getYLims(hParam1TFPlot));
            rescaleData(hParam2TFPlot,xMin,xMax,getYLims(hParam2TFPlot));
            rescaleData(hParam3TFPlot,xMin,xMax,getYLims(hParam3TFPlot));
            rescaleData(hParam4TFPlot,xMin,xMax,getYLims(hParam4TFPlot));
            rescaleData(hParam5TFPlot,xMin,xMax,getYLims(hParam5TFPlot));
            rescaleData(hParam6TFPlot,xMin,xMax,getYLims(hParam6TFPlot));

            for i = 1:size(tfplotHandles,1)
                for j = 1:size(tfplotHandles,2)
                    rescaleData(tfplotHandles(i,j),xMin,xMax,getYLims(tfplotHandles(i,j)));
                end
            end
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType<=3 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        elseif analysisType==6
            xMin = str2double(get(hSTAMin,'String'));
            xMax = str2double(get(hSTAMax,'String'));
        else    
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end

        rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
        rescaleData(hParam1Plot,xMin,xMax,getYLims(hParam1Plot));
        rescaleData(hParam2Plot,xMin,xMax,getYLims(hParam2Plot));
        rescaleData(hParam3Plot,xMin,xMax,getYLims(hParam3Plot));
        rescaleData(hParam4Plot,xMin,xMax,getYLims(hParam4Plot));
        rescaleData(hParam5Plot,xMin,xMax,getYLims(hParam5Plot));
        rescaleData(hParam6Plot,xMin,xMax,getYLims(hParam6Plot));
        
        % [Vinay] - rescale the TF plots
        figure(2);
        
        if plotStyle ~= 1 % pcolor, imagesc
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
            yLims = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
            rescaleData(hParam1TFPlot,xMin,xMax,yLims);
            rescaleData(hParam2TFPlot,xMin,xMax,yLims);
            rescaleData(hParam3TFPlot,xMin,xMax,yLims);
            rescaleData(hParam4TFPlot,xMin,xMax,yLims);
            rescaleData(hParam5TFPlot,xMin,xMax,yLims);
            rescaleData(hParam6TFPlot,xMin,xMax,yLims);

            for i = 1:size(tfplotHandles,1)
                for j = 1:size(tfplotHandles,2)
                    rescaleData(tfplotHandles(i,j),xMin,xMax,yLims);
                end
            end
            
        else % line
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
            rescaleData(hParam1TFPlot,xMin,xMax,getYLims(hParam1TFPlot));
            rescaleData(hParam2TFPlot,xMin,xMax,getYLims(hParam2TFPlot));
            rescaleData(hParam3TFPlot,xMin,xMax,getYLims(hParam3TFPlot));
            rescaleData(hParam4TFPlot,xMin,xMax,getYLims(hParam4TFPlot));
            rescaleData(hParam5TFPlot,xMin,xMax,getYLims(hParam5TFPlot));
            rescaleData(hParam6TFPlot,xMin,xMax,getYLims(hParam6TFPlot));

            for i = 1:size(tfplotHandles,1)
                for j = 1:size(tfplotHandles,2)
                    rescaleData(tfplotHandles(i,j),xMin,xMax,getYLims(tfplotHandles(i,j)));
                end
            end
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');
        
        holdOnGivenPlotHandle(plotHandles,holdOnState);
        holdOnGivenPlotHandle(hParam1Plot,holdOnState);
        holdOnGivenPlotHandle(hParam2Plot,holdOnState);
        holdOnGivenPlotHandle(hParam3Plot,holdOnState);
        holdOnGivenPlotHandle(hParam4Plot,holdOnState);
        holdOnGivenPlotHandle(hParam5Plot,holdOnState);
        holdOnGivenPlotHandle(hParam6Plot,holdOnState);
        
        % TF plots
        holdOnGivenPlotHandle(tfplotHandles,holdOnState);
        holdOnGivenPlotHandle(hParam1TFPlot,holdOnState);
        holdOnGivenPlotHandle(hParam2TFPlot,holdOnState);
        holdOnGivenPlotHandle(hParam3TFPlot,holdOnState);
        holdOnGivenPlotHandle(hParam4TFPlot,holdOnState);
        holdOnGivenPlotHandle(hParam5TFPlot,holdOnState);
        holdOnGivenPlotHandle(hParam6TFPlot,holdOnState);
        
        if holdOnState
            set(hElectrodes,'Nextplot','add');
        else
            set(hElectrodes,'Nextplot','replace');
        end

        function holdOnGivenPlotHandle(plotHandles,holdOnState)
            
            [numRows,numCols] = size(plotHandles);
            if holdOnState
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','add');

                    end
                end
            else
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','replace');
                    end
                end
            end
        end 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        
        claGivenPlotHandle(plotHandles);
        claGivenPlotHandle(hParam1Plot);
        claGivenPlotHandle(hParam2Plot);
        claGivenPlotHandle(hParam3Plot);
        claGivenPlotHandle(hParam4Plot);
        claGivenPlotHandle(hParam5Plot);
        claGivenPlotHandle(hParam6Plot);
        
        cla(hRFMapPlot);cla(hcenterRFMapPlot);
        
        % TF plots
        claGivenPlotHandle(tfplotHandles);
        claGivenPlotHandle(hParam1TFPlot);
        claGivenPlotHandle(hParam2TFPlot);
        claGivenPlotHandle(hParam3TFPlot);
        claGivenPlotHandle(hParam4TFPlot);
        claGivenPlotHandle(hParam5TFPlot);
        claGivenPlotHandle(hParam6TFPlot);
        
        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for i=1:numRows
                for j=1:numCols
                    cla(plotHandles(i,j));
                end
            end
        end
    end
    
    function notch_Callback(hObject,~,~)
       notchData = get(hObject,'val');
    end

    function useDefParams_Callback(hObject,~,~)
       useDefParams = get(hObject,'val');
       if useDefParams == 1
           switch protocolNumber
               case 1 % Ring protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',1); % 1 => NA
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',9); % 9 => rad 
                   set(hParam6,'val',1); % 1 => NA
                   set(hParam7,'val',1); % 1 => NA
               case 2 % Contrast Ring Protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',7); % 7 => cont (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 3 % Dual Contrast Protocol
                   set(hGabor1,'val',2);% 2 => C
                   set(hGabor2,'val',2);
                   set(hGabor3,'val',2);
                   set(hGabor4,'val',2);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',9); % 9 => rad 
                   set(hParam5,'val',7); % 7 => cont (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 4 % Dual Orientation Protocol
                   set(hGabor1,'val',2);% 2 => C
                   set(hGabor2,'val',2);
                   set(hGabor3,'val',2);
                   set(hGabor4,'val',2);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',9); % 9 => rad 
                   set(hParam5,'val',6); % 6 => ori (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 5 % Dual Phase Protocol
                   set(hGabor1,'val',2);% 2 => C
                   set(hGabor2,'val',2);
                   set(hGabor3,'val',2);
                   set(hGabor4,'val',2);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',9); % 9 => rad 
                   set(hParam5,'val',10); % 10 => phase (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 6 %Orientation Ring Protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',6); % 6 => ori (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 7 %Phase Ring Protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',7); % 10 => phase (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 8 %Drifting Ring Protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',9); % 9 => rad
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',8); % 8 => TF (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 9 %Cross Orientation Protocol
                   set(hGabor1,'val',4); % 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',3); % 3 => R
                   set(hGabor5,'val',3);
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',9); % 9 => rad 
                   set(hParam6,'val',9); % 9 => rad (Ring)
                   set(hParam7,'val',1); % 1 => NA
               case 10 %Annulus Fixed Protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',9); % 9 => rad 
                   set(hParam5,'val',7); % 7 => cont (Ring) 
                   set(hParam6,'val',9); % 9 => rad (Ring)
                   set(hParam7,'val',1); % 1 => NA
               otherwise
                   set(hGabor1,'val',4); % 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',3); % 3 => R
                   set(hGabor4,'val',3);
                   set(hGabor5,'val',2); % 2 => C
                   set(hGabor6,'val',2);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',7); % 7 => cont (surr)
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont (ring)
                   set(hParam4,'val',9); % 9 => rad 
                   set(hParam5,'val',7); % 7 => cont (centre)
                   set(hParam6,'val',6); % 6 => ori
                   set(hParam7,'val',1); % 1 => NA
           end
       end
       adjFactor = 5;
       gabor1 = adjFactor - get(hGabor1,'val'); % If sel is C, then value = 2 (NA = 1, C = 2, R = 3, S = 4) => gabor1 = 5 - 2 = 3. If R then 2, if S then 1
       gabor2 = adjFactor - get(hGabor2,'val');
       gabor3 = adjFactor - get(hGabor3,'val');
       gabor4 = adjFactor - get(hGabor4,'val');
       gabor5 = adjFactor - get(hGabor5,'val');
       gabor6 = adjFactor - get(hGabor6,'val');
       gabor7 = adjFactor - get(hGabor7,'val');
       
       param1 = get(hParam1, 'val');
       param2 = get(hParam2, 'val');
       param3 = get(hParam3, 'val');
       param4 = get(hParam4, 'val');
       param5 = get(hParam5, 'val');
       param6 = get(hParam6, 'val');
       param7 = get(hParam7, 'val');
       
       figure(1);
       hParam1Plot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1Grid,0.002);
       hParam2Plot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2Grid,0.002);
       hParam3Plot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3Grid,0.002);
       hParam4Plot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4Grid,0.002);
       hParam5Plot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5Grid,0.002);
       hParam6Plot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6Grid,0.002);
       
       % [Vinay] - reset the TF grid as well
       figure(2);
       hParam1TFPlot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1TFGrid,0.002);
       hParam2TFPlot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2TFGrid,0.002);
       hParam3TFPlot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3TFGrid,0.002);
       hParam4TFPlot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4TFGrid,0.002);
       hParam5TFPlot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5TFGrid,0.002);
       hParam6TFPlot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6TFGrid,0.002);
       
    end

    function resetParams_Callback(~,~)
        adjFactor = 5;
        gabor1 = adjFactor - get(hGabor1,'val'); % If sel is C, then value = 2 (NA = 1, C = 2, R = 3, S = 4) => gabor1 = 5 - 2 = 3. If R then 2, if S then 1
        gabor2 = adjFactor - get(hGabor2,'val');
        gabor3 = adjFactor - get(hGabor3,'val');
        gabor4 = adjFactor - get(hGabor4,'val');
        gabor5 = adjFactor - get(hGabor5,'val');
        gabor6 = adjFactor - get(hGabor6,'val');
        gabor7 = adjFactor - get(hGabor7,'val');
        
        param1 = get(hParam1, 'val');
        param2 = get(hParam2, 'val');
        param3 = get(hParam3, 'val');
        param4 = get(hParam4, 'val');
        param5 = get(hParam5, 'val');
        param6 = get(hParam6, 'val');
        param7 = get(hParam7, 'val');
        
        figure(1);
        hParam1Plot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1Grid,0.002);
        hParam2Plot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2Grid,0.002);
        hParam3Plot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3Grid,0.002);
        hParam4Plot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4Grid,0.002);
        hParam5Plot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5Grid,0.002);
        hParam6Plot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6Grid,0.002);
        
        % [Vinay] - reset the TF grid as well
        figure(2);
        hParam1TFPlot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1TFGrid,0.002);
        hParam2TFPlot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2TFGrid,0.002);
        hParam3TFPlot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3TFGrid,0.002);
        hParam4TFPlot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4TFGrid,0.002);
        hParam5TFPlot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5TFGrid,0.002);
        hParam6TFPlot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6TFGrid,0.002);
        
    end

    function valsUnique = getValsUnique(gaborNumber, paramNumber)
        if gaborNumber > 3
            valsUnique = [];
        else
            switch (paramNumber-1)
                case 1
                    valsUnique = aValsUnique{gaborNumber};
                case 2
                    valsUnique = eValsUnique{gaborNumber};
                case 3
                    valsUnique = sValsUnique{gaborNumber};
                case 4
                    valsUnique = fValsUnique{gaborNumber};
                case 5
                    valsUnique = oValsUnique{gaborNumber};
                case 6
                    valsUnique = cValsUnique{gaborNumber};
                case 7
                    valsUnique = tValsUnique{gaborNumber};
                case 8
                    valsUnique = rValsUnique{gaborNumber};
                case 9
                    valsUnique = pValsUnique{gaborNumber};
                otherwise
                    valsUnique = [];
            end
        end
    end

    function bipolar_Callback(hObject,~,~)
       if (get(hObject,'val') ~= 1)
           useBipolar = 1;
       else
           useBipolar = 0;
       end
    end
        

%---------------TF analysis Functions----------------------------
    function resetMTMParams_Callback(~,~)
        mtmParams.Fs = str2double(get(hMTMFs,'String'));
        mtmParams.tapers(1) = str2double(get(hMTMTapersTW,'String'));
        mtmParams.tapers(2) = str2double(get(hMTMTapersK,'String'));
        
        movingWin(1) = str2double(get(hMTMwLen,'String'));
        movingWin(2) = str2double(get(hMTMwStep,'String'));
    end

    function resetMPParams_Callback(~,~)
        numAtomsMP = str2double(get(hMPnumAtoms,'String'));
    end

    function resetTFSettings_Callback(~,~)
        plotStyle = get(hTFPlotStyle,'val');
        spectrumType = get(hTFSpectrumType,'val');
        
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        tfMethod = get(hTFPlotMethod,'val');
    end

    function rescaleCaxis_Callback(~,~)
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        figure(2);
        set(hParam1TFPlot,'cLim',[cmin cmax]);
        set(hParam2TFPlot,'cLim',[cmin cmax]);
        set(hParam3TFPlot,'cLim',[cmin cmax]);
        set(hParam4TFPlot,'cLim',[cmin cmax]);
        set(hParam5TFPlot,'cLim',[cmin cmax]);
        set(hParam6TFPlot,'cLim',[cmin cmax]);
        
        for i = 1:size(tfplotHandles,1)
            for j = 1:size(tfplotHandles,2)
                set(tfplotHandles(i,j),'cLim',[cmin cmax]);
            end
        end
        
    end
        
        
        

    %----main TF plotting function----------------
    function plotTFData_Callback(~,~)
        a=cell2mat(get(hAzimuth,'val'));
        e=cell2mat(get(hElevation,'val'));
        s=cell2mat(get(hSigma,'val'));
        f=cell2mat(get(hSpatialFreq,'val'));
        o=cell2mat(get(hOrientation,'val'));
        c=cell2mat(get(hContrast,'val'));
        t=cell2mat(get(hTemporalFreq,'val'));
        r=cell2mat(get(hRadius,'val'));
        p=cell2mat(get(hSpatialPhase,'val'));
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        BLMin = str2double(get(hBaselineMin,'String'));
        BLMax = str2double(get(hBaselineMax,'String'));
        STMin = str2double(get(hStimPeriodMin,'String'));
        STMax = str2double(get(hStimPeriodMax,'String'));
        STAMin = str2double(get(hSTAMin,'String'));
        STAMax = str2double(get(hSTAMax,'String'));
        holdOnState = get(hHoldOn,'val');
        removeMeanSTA = get(hRemoveMeanSTA,'val');
        notchData = get(hNotchData, 'val');
        
        
        % Load TF Parameters
        % MTM
        mtmParams.Fs = str2double(get(hMTMFs,'String'));
        mtmParams.tapers(1) = str2double(get(hMTMTapersTW,'String'));
        mtmParams.tapers(2) = str2double(get(hMTMTapersK,'String'));
        movingWin(1) = str2double(get(hMTMwLen,'String'));
        movingWin(2) = str2double(get(hMTMwStep,'String'));
        
        % MP
        numAtomsMP = str2double(get(hMPnumAtoms,'String'));
        
        % Plot settings
        plotStyle = get(hTFPlotStyle,'val');
        spectrumType = get(hTFSpectrumType,'val');
        
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        tfMethod = get(hTFPlotMethod,'val');
        
        
        % [Vinay] - get the plot handles here again
        figure(2);
        hParam1TFPlot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1TFGrid,0.002);
        hParam2TFPlot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2TFGrid,0.002);
        hParam3TFPlot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3TFGrid,0.002);
        hParam4TFPlot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4TFGrid,0.002);
        hParam5TFPlot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5TFGrid,0.002);
        hParam6TFPlot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6TFGrid,0.002);
        
        
        analogChannelPos = get(hAnalogChannel,'val');
        analogChannelString = analogChannelStringArray{analogChannelPos};

        % [Vinay] - read the analog channel 2 string for bipolar case
        analogChannelPos2 = get(hAnalogChannel2,'val');
        analogChannelString2 = ('none');
        if analogChannelPos2 ~= 1
            analogChannelString2 = analogChannelStringArray{analogChannelPos2-1};
        end

        tfplotLFPData1Channel(tfplotHandles,analogChannelString,analogChannelString2,s,f,o,c,t,r,p,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState);
        tfplotLFPData1Parameter1Channel(hParam1TFPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor1,param1,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState);
        tfplotLFPData1Parameter1Channel(hParam2TFPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor2,param2,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState);
        tfplotLFPData1Parameter1Channel(hParam3TFPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor3,param3,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState);
        tfplotLFPData1Parameter1Channel(hParam4TFPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor4,param4,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState);
        tfplotLFPData1Parameter1Channel(hParam5TFPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor5,param5,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState);
        tfplotLFPData1Parameter1Channel(hParam6TFPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor6,param6,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState);
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that plots the data
function rfMapVals = plotLFPData1Channel(plotHandles,channelString,analogChannelString2,s,f,o,c,t,r,p,folderLFP,...
analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName,notchData,useBipolar) % Vinay - added notchData

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

titleFontSize = 10;

[parameterCombinations,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

% Get the data
removeAvgRef = 0;
if removeAvgRef
    disp('Removing average reference');
    load(fullfile(folderLFP,'avgRef'));
    avgRef = analogData;
end
clear signal analogData
load(fullfile(folderLFP,channelString));
if removeAvgRef
    analogData = analogData-avgRef;
end

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

rfMapVals = zeros(numRows,numCols);
for i=1:numRows
    e = numRows-i+1;
    for j=1:numCols
        a = j;
        clear goodPos
%         goodPos = parameterCombinations{a,e,s,f,o,c,t,r,p};
        % [Vinay] - goodPos will be obtained by taking the goodPos for all
        % the three gabors and finding the common trials in them
        clear pos1 pos2 pos3
        pos1 = parameterCombinations{a,e,s(1),f(1),o(1),c(1),t(1),p(1),r(1),1}; % good positions for gabor 1 i.e. S
        pos2 = parameterCombinations{a,e,s(2),f(2),o(2),c(2),t(2),p(2),r(2),2}; % for R
        pos3 = parameterCombinations{a,e,s(3),f(3),o(3),c(3),t(3),p(3),r(3),3}; % for C
        goodPos = intersect(pos1,pos2);
        goodPos = intersect(goodPos,pos3);
        goodPos = setdiff(goodPos,badTrials);
      
        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            disp(['pos=(' num2str(i) ',' num2str(j) ') ,n=' num2str(length(goodPos))]);
    
            Fs = round(1/(timeVals(2)-timeVals(1)));
            BLRange = uint16((BLMax-BLMin)*Fs);
            STRange = uint16((STMax-STMin)*Fs);
            BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
            STPos = find(timeVals>=STMin,1)+ (1:STRange);

            xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
            xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);
            
%             % Vinay - added next two 'if' conditions to take care of the case
%             % when the lengths of BLPos(STPos) and xsBL(xsST) differ by 1
%             if ((abs(length(BLRange)-length(xsBL)))==1)
%                 BLRange = ceil((BLMax-BLMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
%             end
%         
%             if ((abs(length(STRange)-length(xsST)))==1)
%                 STRange = ceil((STMax-STMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
%             end
            
            if notchData
                analogData = analogDataNotched;
            end
            
            if useBipolar
                analogChannelString1 = channelString; % first electrode selected
                clear analogData analogDataNotched
                load(fullfile(folderLFP,analogChannelString1));
                
                if notchData
                    analogData1 = analogDataNotched;
                else
                    analogData1 = analogData; % Vinay - these are the 
                    % freshly loaded values for electrode 1
                end
                % Vinay - store the analogData 
                % for electrode 1 in a separate variable 
                % otherwise loading the data for 2nd electrode will
                % overwrite analogData
                
                clear analogData analogDataNotched
                load(fullfile(folderLFP,analogChannelString2));
                
                if notchData
                    analogData2 = analogDataNotched;
                else
                    analogData2 = analogData; % Vinay - these are the 
                    % freshly loaded values for electrode 2
                end
                analogData = analogData1 - analogData2;
            end
                

            if analysisType == 1        % compute ERP
                clear erp
                erp = mean(analogData(goodPos,:),1); %#ok<*NODEF>
                plot(plotHandles(i,j),timeVals,erp,'color',plotColor);
                
                rfMapVals(e,a) = rms(erp(STPos));

            elseif analysisType == 2  ||   analysisType == 3 % compute Firing rates
                disp('Use plotSpikeData instead of plotLFPData...');
            else
                
                fftBL = abs(fft(analogData(goodPos,BLPos),[],2));
                fftST = abs(fft(analogData(goodPos,STPos),[],2));

                if analysisType == 4
                    plot(plotHandles(i,j),xsBL,log10(mean(fftBL)),'g');
                    set(plotHandles(i,j),'Nextplot','add');
                    plot(plotHandles(i,j),xsST,log10(mean(fftST)),'k');
                    set(plotHandles(i,j),'Nextplot','replace');
                end

                if analysisType == 5
                    if xsBL == xsST %#ok<BDSCI>
                        plot(plotHandles(i,j),xsBL,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor);
                    else
                        disp('Choose same baseline and stimulus periods..');
                    end
                end
            end
            
            % Display title
            if (i==1)
                if (j==1)
                    title(plotHandles(i,j),['Azi: ' num2str(aValsUnique{1}(a))],'FontSize',titleFontSize);
                    % [Vinay] - because all the gabors are concentric,
                    % using the aValsUnique for the gabor 1 w.l.g.
                else
                    title(plotHandles(i,j),num2str(aValsUnique{1}(a)),'FontSize',titleFontSize);
                end
            end
                
            if (j==numCols)
                if (i==1)
                 title(plotHandles(i,j),[{'Ele'} {num2str(eValsUnique{1}(e))}],'FontSize',titleFontSize,...
                     'Units','Normalized','Position',[1.25 0.5]);
                 % [Vinay] - because all the gabors are concentric,
                    % using the eValsUnique for the gabor 1 w.l.g.
                else
                    title(plotHandles(i,j),num2str(eValsUnique{1}(e)),'FontSize',titleFontSize,...
                     'Units','Normalized','Position',[1.25 0.5]);
                end
            end
        end
    end
end

if analysisType~=1
    rfMapVals=[];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLFPData1Parameter1Channel(plotHandles,channelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gaborNum,paramNum,folderLFP,...
analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar) % Vinay - added notchData

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');


titleFontSize = 10;

timeForComputation = [40 100]/1000; % ms
freqForComputation = [40 60]; % Hz

[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);
[~,numCols] = size(plotHandles);

% Get the data
clear signal analogData
load(fullfile(folderLFP,channelString));

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

% [Vinay] - repeat the set of parameters for each gabor depending on the
% number of columns to be drawn
aList = repmat(a,1,numCols);
eList = repmat(e,1,numCols);
sList = repmat(s,1,numCols);
fList = repmat(f,1,numCols);
oList = repmat(o,1,numCols);
cList = repmat(c,1,numCols);
tList = repmat(t,1,numCols);
rList = repmat(r,1,numCols);
pList = repmat(p,1,numCols);

% [Vinay] decide the parameter based on paramNum
switch (paramNum-1)
    case 1
        aList(gaborNum,:) = 1:length(aValsUnique{gaborNum}); 
        titleParam = 'Azi: ';
        titleList = aValsUnique{gaborNum};
    case 2
        eList(gaborNum,:) = 1:length(eValsUnique{gaborNum}); 
        titleParam = 'Ele: ';
        titleList = eValsUnique{gaborNum};
    case 3
        sList(gaborNum,:) = 1:length(sValsUnique{gaborNum}); 
        titleParam = 'Sigma: ';
        titleList = sValsUnique{gaborNum};
    case 4
        fList(gaborNum,:) = 1:length(fValsUnique{gaborNum}); 
        titleParam = 'SF: ';
        titleList = fValsUnique{gaborNum};
    case 5
        oList(gaborNum,:) = 1:length(oValsUnique{gaborNum}); 
        titleParam = 'Ori: ';
        titleList = oValsUnique{gaborNum};
    case 6
        cList(gaborNum,:) = 1:length(cValsUnique{gaborNum}); 
        titleParam = 'Contr: ';
        titleList = cValsUnique{gaborNum};
    case 7
        tList(gaborNum,:) = 1:length(tValsUnique{gaborNum}); 
        titleParam = 'TF: ';
        titleList = tValsUnique{gaborNum};
    case 8
        rList(gaborNum,:) = 1:length(rValsUnique{gaborNum}); 
        titleParam = 'Rad: ';
        titleList = rValsUnique{gaborNum};
    case 9
        pList(gaborNum,:) = 1:length(pValsUnique{gaborNum}); 
        titleParam = 'SPhs: ';
        titleList = pValsUnique{gaborNum};
    otherwise
        disp('No particular gabor or parameter selected');
end

% [Vinay] - Get the lengths of indices in parameterCombinations
for i = 1:3
    aLen(i) = length(aValsUnique{i});
    eLen(i) = length(eValsUnique{i});
    sLen(i) = length(sValsUnique{i});
    fLen(i) = length(fValsUnique{i});
    oLen(i) = length(oValsUnique{i});
    cLen(i) = length(cValsUnique{i});
    tLen(i) = length(tValsUnique{i});
    pLen(i) = length(pValsUnique{i}); % Vinay - for CRS
    rLen(i) = length(rValsUnique{i}); % Vinay - for CRS
    
    % If more than one value, then length is one greater for all the values
    % together
    if (aLen(i)> 1)           aLen(i)=aLen(i)+1;                    end
    if (eLen(i)> 1)           eLen(i)=eLen(i)+1;                    end
    if (sLen(i)> 1)           sLen(i)=sLen(i)+1;                    end
    if (fLen(i)> 1)           fLen(i)=fLen(i)+1;                    end
    if (oLen(i)> 1)           oLen(i)=oLen(i)+1;                    end
    if (cLen(i)> 1)           cLen(i)=cLen(i)+1;                    end
    if (tLen(i)> 1)           tLen(i)=tLen(i)+1;                    end
    if (pLen(i)> 1)           pLen(i)=pLen(i)+1;                    end % Vinay - for CRS
    if (rLen(i)> 1)           rLen(i)=rLen(i)+1;                    end % Vinay - for CRS
end
    
% Main loop
computationVals=zeros(1,numCols);
for j=1:numCols
    clear goodPos
%     goodPos = parameterCombinations{aList(j),eList(j),sList(j),fList(j),oList(j),cList(j),tList(j)};
%     goodPos = setdiff(goodPos,badTrials);
    
    % [Vinay] - goodPos will be obtained by taking the goodPos for all
    % the three gabors and finding the common trials in them
    clear pos1 pos2 pos3
    switch protocolNumber
               case 1 % Ring protocol
                   pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                   pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad is imp
                   % [Vinay] -  C is matched to S and therefore the entries
                   % for gabor3 i.e. can be ignored. Otherwise this can 
                   % give spurious results if the param values for C are 
                   % not set equal to those for S. Hence pos3 is the full
                   % set here
                   pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C
                   
               case 2 % Contrast Ring Protocol
                   pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                   pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2,j),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & cont
                   % [Vinay] -  C is matched to S and therefore the entries
                   % for gabor3 i.e. can be ignored. Otherwise this can 
                   % give spurious results if the param values for C are 
                   % not set equal to those for S. Hence pos3 is the full
                   % set here
                   pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C

               case 3 % Dual Contrast Protocol
                   pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                   % [Vinay] - S is hidden in this case and only C & R
                   % define the stimuli. pos1 is therefore the full set
                   % here
                   pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2,j),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & cont
                   % [Vinay] -  R params - cont & rad. Rest are matched to
                   % C's params
                   pos3 = parameterCombinations{aList(3,j),eList(3,j),sList(3,j),fList(3,j),oList(3,j),cList(3,j),tList(3,j),pList(3,j),rList(3,j),3}; % for C

               case 4 % Dual Orientation Protocol
                   pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                   % [Vinay] - S is hidden in this case and only C & R
                   % define the stimuli. pos1 is therefore the full set
                   % here
                   pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2,j),cLen(2),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & ori
                   % [Vinay] -  R params - ori & rad. Rest are matched to
                   % C's params
                   pos3 = parameterCombinations{aList(3,j),eList(3,j),sList(3,j),fList(3,j),oList(3,j),cList(3,j),tList(3,j),pList(3,j),rList(3,j),3}; % for C

               case 5 % Dual Phase Protocol
                   pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                   % [Vinay] - S is hidden in this case and only C & R
                   % define the stimuli. pos1 is therefore the full set
                   % here
                   pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2,j),rList(2,j),2}; % for R - only rad & p
                   % [Vinay] -  R params - p & rad. Rest are matched to
                   % C's params
                   pos3 = parameterCombinations{aList(3,j),eList(3,j),sList(3,j),fList(3,j),oList(3,j),cList(3,j),tList(3,j),pList(3,j),rList(3,j),3}; % for C

               case 6 %Orientation Ring Protocol
                   pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                   pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2,j),cLen(2),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & ori
                   % [Vinay] -  C is matched to S and therefore the entries
                   % for gabor3 i.e. can be ignored. Otherwise this can 
                   % give spurious results if the param values for C are 
                   % not set equal to those for S. Hence pos3 is the full
                   % set here
                   pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C

               case 7 %Phase Ring Protocol
                   pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                   pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2,j),rList(2,j),2}; % for R - only rad & p
                   % [Vinay] -  C is matched to S and therefore the entries
                   % for gabor3 i.e. can be ignored. Otherwise this can 
                   % give spurious results if the param values for C are 
                   % not set equal to those for S. Hence pos3 is the full
                   % set here
                   pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C

               case 8 %Drifting Ring Protocol
                   pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                   pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tList(2,j),pLen(2),rList(2,j),2}; % for R - only rad & t
                   % [Vinay] -  C is matched to S and therefore the entries
                   % for gabor3 i.e. can be ignored. Otherwise this can 
                   % give spurious results if the param values for C are 
                   % not set equal to those for S. Hence pos3 is the full
                   % set here
                   pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C

               case 9 %Cross Orientation Protocol
                   disp('Left for future extension :)');
               case 10 %Annulus Fixed Protocol
                   pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                   pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2,j),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & cont
                   % [Vinay] -  C is matched to S and therefore the entries
                   % for gabor3 i.e. can be ignored. Otherwise this can 
                   % give spurious results if the param values for C are 
                   % not set equal to those for S. However, the rad of C is
                   % still a relevant parameter and has to be taken care of
                   pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3,j),3}; % for C - only rad
                   
        otherwise
            pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1};
            % good positions for gabor 1 i.e. S
            pos2 = parameterCombinations{aList(2,j),eList(2,j),sList(2,j),fList(2,j),oList(2,j),cList(2,j),tList(2,j),pList(2,j),rList(2,j),2}; % for R
            pos3 = parameterCombinations{aList(3,j),eList(3,j),sList(3,j),fList(3,j),oList(3,j),cList(3,j),tList(3,j),pList(3,j),rList(3,j),3}; % for C
    end
    
    goodPos = intersect(pos1,pos2);
    goodPos = intersect(goodPos,pos3);
    goodPos = setdiff(goodPos,badTrials);

    if isempty(goodPos)
        disp('No entries for this combination..')
    else
        disp(['pos=' num2str(j) ',n=' num2str(length(goodPos))]);

        Fs = round(1/(timeVals(2)-timeVals(1)));
        BLRange = uint16((BLMax-BLMin)*Fs);
        STRange = uint16((STMax-STMin)*Fs);
        BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
        STPos = find(timeVals>=STMin,1)+ (1:STRange);

        xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
        xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);
        
%         % Vinay - added next two 'if' conditions to take care of the case
%         % when the lengths of BLPos(STPos) and xsBL(xsST) differ by 1
%         if ((abs(length(BLRange)-length(xsBL)))==1)
%             BLRange = ceil((BLMax-BLMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
%         end
%         
%         if ((abs(length(STRange)-length(xsST)))==1)
%             STRange = ceil((STMax-STMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
%         end
        
        xsComputation = intersect(find(timeVals>=timeForComputation(1)),find(timeVals<timeForComputation(2)));
        freqComputation = intersect(find(xsST>=freqForComputation(1)),find(xsST<=freqForComputation(2)));
        
        % Vinay - added this notch data check
        if notchData
            analogData = analogDataNotched;
        end
        
        if useBipolar
                analogChannelString1 = channelString; % first electrode selected
                clear analogData analogDataNotched
                load(fullfile(folderLFP,analogChannelString1));
                
                if notchData
                    analogData1 = analogDataNotched;
                else
                    analogData1 = analogData; % Vinay - these are the 
                    % freshly loaded values for electrode 1
                end
                % Vinay - store the analogData 
                % for electrode 1 in a separate variable 
                % otherwise loading the data for 2nd electrode will
                % overwrite analogData
                
                clear analogData analogDataNotched
                load(fullfile(folderLFP,analogChannelString2));
                
                if notchData
                    analogData2 = analogDataNotched;
                else
                    analogData2 = analogData; % Vinay - these are the 
                    % freshly loaded values for electrode 2
                end
                analogData = analogData1 - analogData2;
        end
        
        if analysisType == 1        % compute ERP
            clear erp
            erp = mean(analogData(goodPos,:),1);
            plot(plotHandles(j),timeVals,erp,'color',plotColor);
            
%             if isempty(o) % Orientation tuning
%                 computationVals(j) = abs(min(erp(xsComputation)));
%             end
            
            if paramNum==6 % Orientation tuning
                computationVals(j) = abs(min(erp(xsComputation)));
            end


        elseif analysisType == 2 || analysisType == 3   % compute Firing rates
            disp('Use plotSpikeData instead of plotLFPData...');
        else

            fftBL = abs(fft(analogData(goodPos,BLPos),[],2));
            fftST = abs(fft(analogData(goodPos,STPos),[],2));

            if analysisType == 4
                plot(plotHandles(j),xsBL,log10(mean(fftBL)),'g');
                set(plotHandles(j),'Nextplot','add');
                plot(plotHandles(j),xsST,log10(mean(fftST)),'k');
                set(plotHandles(j),'Nextplot','replace');
            end

            if analysisType == 5
                if xsBL == xsST %#ok<BDSCI>
                    plot(plotHandles(j),xsBL,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor);
                else
                    disp('Choose same baseline and stimulus periods..');
                end
            end
            
%             if isempty(o) % Orientation tuning
%                 computationVals(j) = max(mean(fftST(:,freqComputation),1));
%             end
            if paramNum==6 % Orientation tuning
                computationVals(j) = max(mean(fftST(:,freqComputation),1));
            end
        end

        % Display title
        if (j==1)
            title(plotHandles(j),[titleParam num2str(titleList(j))],'FontSize',titleFontSize);
        else
            title(plotHandles(j),num2str(titleList(j)),'FontSize',titleFontSize);
        end
    end
end

% Orientation tuning
if paramNum==6
    disp(['o: ' num2str(computationVals)]);
    [prefOrientation,orientationSelectivity] = getOrientationTuning(computationVals,oValsUnique{gaborNum});
    disp(['prefOri: ' num2str(round(prefOrientation)) ', sel: ' num2str(orientationSelectivity)]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSpikeData1Channel(plotHandles,channelNumber,s,f,o,c,t,folderSpikes,...
analysisType,timeVals,plotColor,unitID,folderName)
titleFontSize = 12;

if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

[parameterCombinations,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

% Get the data
clear signal spikeData
load([folderSpikes 'elec' num2str(channelNumber) '_SID' num2str(unitID)]);

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

for i=1:numRows
    e = numRows-i+1;
    for j=1:numCols
        a = j;
        clear goodPos
        goodPos = parameterCombinations{a,e,s,f,o,c,t};
        goodPos = setdiff(goodPos,badTrials);

        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            disp(['pos=(' num2str(i) ',' num2str(j) ') ,n=' num2str(length(goodPos))]);
            
            if analysisType == 2
                [psthVals,xs] = getPSTH(spikeData(goodPos),10,timeVals(1),timeVals(end));
                plot(plotHandles(i,j),xs,psthVals,'color',plotColor);
            else
                X = spikeData(goodPos);
                axes(plotHandles(i,j)); %#ok<LAXES>
                rasterplot(X,1:length(X),plotColor);
            end
        end
        
        % Display title
        if (i==1)
            if (j==1)
                title(plotHandles(i,j),['Azi: ' num2str(aValsUnique(a))],'FontSize',titleFontSize);
            else
                title(plotHandles(i,j),num2str(aValsUnique(a)),'FontSize',titleFontSize);
            end
        end

        if (j==numCols)
            if (i==1)
                title(plotHandles(i,j),[{'Ele'} {num2str(eValsUnique(e))}],'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            else
                title(plotHandles(i,j),num2str(eValsUnique(e)),'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSpikeData1Parameter1Channel(plotHandles,channelNumber,a,e,s,f,o,c,t,folderSpikes,...
analysisType,timeVals,plotColor,unitID,folderName)
titleFontSize = 12;

if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);
[~,numCols] = size(plotHandles);

% Get the data
clear signal spikeData
load([folderSpikes 'elec' num2str(channelNumber) '_SID' num2str(unitID)]);

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

% Out of a,e,s,f,o,c ant t, only one parameter is empty
if isempty(a)
    aList = 1:length(aValsUnique); 
    eList = e+zeros(1,numCols); 
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Azi: ';
    titleList = aValsUnique;
end

if isempty(e) 
    aList = a+zeros(1,numCols);
    eList = 1:length(eValsUnique);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Ele: ';
    titleList = eValsUnique;
end

if isempty(s) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = 1:length(sValsUnique);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Sigma: ';
    titleList = sValsUnique;
end

if isempty(f)
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = 1:length(fValsUnique); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'SF: ';
    titleList = fValsUnique;
end

if isempty(o) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = 1:length(oValsUnique);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Ori: ';
    titleList = oValsUnique;
end

if isempty(c) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = 1:length(cValsUnique);
    tList = t+zeros(1,numCols);
    titleParam = 'Con: ';
    titleList = cValsUnique;
end

if isempty(t) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols);
    tList = 1:length(tValsUnique);
    titleParam = 'TF: ';
    titleList = tValsUnique;
end

% Plot

for j=1:numCols
    %a = j;
    clear goodPos
    goodPos = parameterCombinations{aList(j),eList(j),sList(j),fList(j),oList(j),cList(j),tList(j)};
    goodPos = setdiff(goodPos,badTrials);

    if isempty(goodPos)
        disp('No entries for this combination..')
    else
        disp(['pos=' num2str(j) ',n=' num2str(length(goodPos))]);
        if analysisType == 2
            [psthVals,xs] = getPSTH(spikeData(goodPos),10,timeVals(1),timeVals(end));
            plot(plotHandles(j),xs,psthVals,'color',plotColor);
        else
            X = spikeData(goodPos);
            axes(plotHandles(j)); %#ok<LAXES>
            rasterplot(X,1:length(X),plotColor);
        end
    end

    % Display title
    if (j==1)
        title(plotHandles(j),[titleParam num2str(titleList(j))],'FontSize',titleFontSize);
    else
        title(plotHandles(j),num2str(titleList(j)),'FontSize',titleFontSize);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSTA1Channel(plotHandles,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
                s,f,o,c,t,timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName,staLen,removeMeanSTA)

titleFontSize = 12;

if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

[parameterCombinations,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

% Get the analog data
clear signal analogData
load([folderLFP analogChannelString]);

% Get the spike data
clear signal spikeData
load([folderSpikes 'elec' num2str(spikeChannelNumber) '_SID' num2str(unitID)]);

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

staTimeLims{1} = [BLMin BLMax];
staTimeLims{2} = [STMin STMax];

for i=1:numRows
    e = numRows-i+1;
    for j=1:numCols
        a = j;
        clear goodPos
        goodPos = parameterCombinations{a,e,s,f,o,c,t};
        goodPos = setdiff(goodPos,badTrials);

        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            goodSpikeData = spikeData(goodPos);
            goodAnalogSignal = analogData(goodPos,:);
            [staVals,numberOfSpikes,xsSTA] = getSTA(goodSpikeData,goodAnalogSignal,staTimeLims,timeVals,staLen,removeMeanSTA);
            
            disp([i j ', numStim: ' length(goodPos) ', numSpikes: ' num2str(numberOfSpikes)]);
            plot(plotHandles(i,j),xsSTA,staVals{1},'color',plotColors{1});
            set(plotHandles(i,j),'Nextplot','add');
            plot(plotHandles(i,j),xsSTA,staVals{2},'color',plotColors{2});
            set(plotHandles(i,j),'Nextplot','replace');
        end
        
        % Display title
        if (i==1)
            if (j==1)
                title(plotHandles(i,j),['Azi: ' num2str(aValsUnique(a))],'FontSize',titleFontSize);
            else
                title(plotHandles(i,j),num2str(aValsUnique(a)),'FontSize',titleFontSize);
            end
        end

        if (j==numCols)
            if (i==1)
                title(plotHandles(i,j),[{'Ele'} {num2str(eValsUnique(e))}],'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            else
                title(plotHandles(i,j),num2str(eValsUnique(e)),'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plotSTA1Parameter1Channel(hParam3Plot,analogChannelNumber,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
%                 a,e,s,f,o,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotRFMaps(hRFMapPlot,hcenterRFMapPlot,rfMapVals,aValsUnique,eValsUnique,plotColor,holdOnState)

plotSimple = 0;% Just compute the mean and Variance

if plotSimple
    [aziCenter,eleCenter] = getRFcenterSimple(aValsUnique,eValsUnique,rfMapVals); %#ok<*UNRCH>
else
    outParams = getRFcenter(aValsUnique,eValsUnique,rfMapVals);
    aziCenter = outParams(1); eleCenter = outParams(2);
    RFSize = sqrt((outParams(3)^2+outParams(4)^2)/2);
end

if plotSimple
    set(hRFMapPlot,'visible','off')
else
    % Plot the gaussian
    dX = (aValsUnique(end)-aValsUnique(1))/100;
    dY = (eValsUnique(end)-eValsUnique(1))/100;
    [~,outVals,boundaryX,boundaryY] = gauss2D(outParams,aValsUnique(1):dX:aValsUnique(end),eValsUnique(1):dY:eValsUnique(end));
    pcolor(hRFMapPlot,aValsUnique(1):dX:aValsUnique(end),eValsUnique(1):dY:eValsUnique(end),outVals);
    shading(hRFMapPlot,'interp'); 
    set(hRFMapPlot,'Nextplot','add');
    plot(hRFMapPlot,boundaryX,boundaryY,'k');
    set(hRFMapPlot,'Nextplot','replace');
end

% Plot the center only
if holdOnState
    set(hcenterRFMapPlot,'Nextplot','add');
else
    set(hcenterRFMapPlot,'Nextplot','replace');
end
plot(hcenterRFMapPlot,aziCenter,eleCenter,[plotColor '+']);
title(hcenterRFMapPlot,['Loc: (' num2str(aziCenter) ',' num2str(eleCenter) '), Size: ' num2str(RFSize)]);
axis(hcenterRFMapPlot,[aValsUnique(1) aValsUnique(end) eValsUnique(1) eValsUnique(end)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yLims = getYLims(plotHandles)

[numRows,numCols] = size(plotHandles);
% Initialize
yMin = inf;
yMax = -inf;

for row=1:numRows
    for column=1:numCols
        % get positions
        axis(plotHandles(row,column),'tight');
        tmpAxisVals = axis(plotHandles(row,column));
        if tmpAxisVals(3) < yMin
            yMin = tmpAxisVals(3);
        end
        if tmpAxisVals(4) > yMax
            yMax = tmpAxisVals(4);
        end
    end
end

yLims=[yMin yMax];
end
function rescaleData(plotHandles,xMin,xMax,yLims)

[numRows,numCols] = size(plotHandles);
labelSize=12;
for i=1:numRows
    for j=1:numCols
        axis(plotHandles(i,j),[xMin xMax yLims]);
        if (i==numRows && rem(j,2)==1)
            if j~=1
                set(plotHandles(i,j),'YTickLabel',[],'fontSize',labelSize);
            end
        elseif (rem(i,2)==0 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
        else
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        end
    end
end

% Remove Labels on the four corners
%set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outString = getStringFromValues(valsUnique,decimationFactor, gaborNum)

if length(valsUnique{gaborNum})==1
    outString = convertNumToStr(valsUnique{gaborNum}(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique{gaborNum})
        outString = cat(2,outString,[convertNumToStr(valsUnique{gaborNum}(i),decimationFactor) '|']);
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
function [outString,outArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums)
outString='';
count=1;
for i=1:length(analogChannelsStored)
    outArray{count} = ['elec' num2str(analogChannelsStored(i))]; %#ok<AGROW>
    outString = cat(2,outString,[outArray{count} '|']);
    count=count+1;
end
if ~isempty(analogInputNums)
    for i=1:length(analogInputNums)
        outArray{count} = ['ainp' num2str(analogInputNums(i))]; %#ok<AGROW>
        outString = cat(2,outString,[outArray{count} '|']);
        count=count+1;
    end
end
end
function outString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs)
outString='';
for i=1:length(neuralChannelsStored)
    outString = cat(2,outString,[num2str(neuralChannelsStored(i)) ', SID ' num2str(SourceUnitIDs(i)) '|']);
end 
end
function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%
%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Data
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo'));
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end
function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
fileName = [folderSpikes 'spikeInfo.mat'];
if exist(fileName,'file')
    load(fileName);
else
    neuralChannelsStored=[];
    SourceUnitID=[];
end
end
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique, rValsUnique, pValsUnique] = loadParameterCombinations(folderExtract)
% [Vinay] - added rValsUnique and pValsUnique for radius and spatial phase
load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('sValsUnique','var')
    sValsUnique=rValsUnique;
end

if ~exist('cValsUnique','var')
    cValsUnique=[];
end

if ~exist('tValsUnique','var')
    tValsUnique=[];
end
end
% function stimResults = loadStimResults(folderExtract)
% load ([folderExtract 'stimResults']);
% end
function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [aziCenter,eleCenter] = getRFcenterSimple(aValsUnique,eValsUnique,rfMapValues)

rfMapValues = rfMapValues-mean(mean(rfMapValues));
% Compute the mean and variance
[maxEnonUnique,maxAnonUnique] = (find(rfMapValues==max(max(rfMapValues))));
maxE = maxEnonUnique(ceil(length(maxEnonUnique)/2));
maxA = maxAnonUnique(ceil(length(maxAnonUnique)/2));

%[maxE maxA]
aziCenter = sum(rfMapValues(maxE,:).*aValsUnique)/sum(rfMapValues(maxE,:));
eleCenter = sum(rfMapValues(:,maxA)'.*eValsUnique)/sum(rfMapValues(:,maxA));
end
function [prefOrientation,orientationSelectivity] = getOrientationTuning(computationVals,oValsUnique)
num=0;
den=0;

for j=1:length(oValsUnique)
    num = num+computationVals(j)*sind(2*oValsUnique(j));
    den = den+computationVals(j)*cosd(2*oValsUnique(j));
end

prefOrientation = 90*atan2(num,den)/pi;
orientationSelectivity = abs(den+1i*num)/sum(computationVals);

if prefOrientation<0
    prefOrientation = prefOrientation+180;
end
end

function protocolNumber = getProtocolNumber(folderExtract)
load (fullfile(folderExtract,'stimResults'));
protocolNumber = stimResults.protocolNumber;
end

function gaborsDisplayed = getGaborsDisplayed(folderExtract)
load (fullfile(folderExtract,'stimResults'));
gaborsDisplayed = stimResults.side;
end


%--------- TF plots functions

function tfplotLFPData1Channel(tfplotHandles,channelString,analogChannelString2,s,f,o,c,t,r,p,folderLFP,...
timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName,notchData,useBipolar,...
tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax,holdOnState) % Vinay - tf parameters

if ispc
    folderExtract = fullfile(folderName,'extractedData');
    folderSegment = fullfile(folderName,'segmentedData');
    mpFolder = fullfile(folderName,'mpAnalysis'); % for using stored MP data
    mpFolderTemp = 'C:\Documents\MP\data\'; % for online calculations
else
    folderExtract = fullfile(folderName,'extractedData');
    folderSegment = fullfile(folderName,'segmentedData');
    mpFolder = fullfile(folderName,'mpAnalysis'); % for using stored MP data
    mpFolderTemp = '/home/vinay/Documents/MP/Data/'; % for online calculations
end

titleFontSize = 9;

[parameterCombinations,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(tfplotHandles);

% Get the data
removeAvgRef = 0;
if removeAvgRef
    disp('Removing average reference');
    load(fullfile(folderLFP,'avgRef'));
    avgRef = analogData;
end
clear signal analogData
load(fullfile(folderLFP,channelString));
if removeAvgRef
    analogData = analogData-avgRef;
end

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

    for i=1:numRows
        e = numRows-i+1;
        for j=1:numCols
            a = j;
            clear goodPos
    %         goodPos = parameterCombinations{a,e,s,f,o,c,t,r,p};
            % [Vinay] - goodPos will be obtained by taking the goodPos for all
            % the three gabors and finding the common trials in them
            clear pos1 pos2 pos3
            pos1 = parameterCombinations{a,e,s(1),f(1),o(1),c(1),t(1),p(1),r(1),1}; % good positions for gabor 1 i.e. S
            pos2 = parameterCombinations{a,e,s(2),f(2),o(2),c(2),t(2),p(2),r(2),2}; % for R
            pos3 = parameterCombinations{a,e,s(3),f(3),o(3),c(3),t(3),p(3),r(3),3}; % for C
            goodPos = intersect(pos1,pos2);
            goodPos = intersect(goodPos,pos3);
            goodPos = setdiff(goodPos,badTrials);

            if isempty(goodPos)
                disp('No entries for this combination..')
            else
                disp(['pos=(' num2str(i) ',' num2str(j) ') ,n=' num2str(length(goodPos))]);

                Fs = round(1/(timeVals(2)-timeVals(1)));
                BLRange = uint16((BLMax-BLMin)*Fs);
                STRange = uint16((STMax-STMin)*Fs);
                BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
                STPos = find(timeVals>=STMin,1)+ (1:STRange);

                xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
                xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);


                if notchData
                    analogData = analogDataNotched;
                end

                if useBipolar
                    analogChannelString1 = channelString; % first electrode selected
                    clear analogData analogDataNotched
                    load(fullfile(folderLFP,analogChannelString1));

                    if notchData
                        analogData1 = analogDataNotched;
                    else
                        analogData1 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 1
                    end
                    % Vinay - store the analogData 
                    % for electrode 1 in a separate variable 
                    % otherwise loading the data for 2nd electrode will
                    % overwrite analogData

                    clear analogData analogDataNotched
                    load(fullfile(folderLFP,analogChannelString2));

                    if notchData
                        analogData2 = analogDataNotched;
                    else
                        analogData2 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 2
                    end
                    analogData = analogData1 - analogData2;
                end

                %--------------------------------------------------------------
                % Implement the method set by tfMethod varibale

                if (tfMethod == 1) % MTM
                    %----------------------------------------------------------
                    if (plotStyle == 3) % line plot (power vs f)
                        if (spectrumType == 1) % raw
                            [SBL,f]=mtspectrumc((analogData(goodPos,BLPos))',mtmParams); % baseline period
                            [SST,f]=mtspectrumc((analogData(goodPos,STPos))',mtmParams); % stimulus period
                            plot(tfplotHandles(i,j),f,log10(SST),'color',plotColor);
                            set(tfplotHandles(i,j),'Nextplot','add');
                            plot(tfplotHandles(i,j),f,log10(SBL),'color','g');
                            if holdOnState
                                set(tfplotHandles(i,j),'Nextplot','add');
                            else
                                set(tfplotHandles(i,j),'Nextplot','replace');
                            end

                        elseif (spectrumType == 2) % difference from baseline
                            if xsST == xsBL 
                                [SBL,f]=mtspectrumc((analogData(goodPos,BLPos))',mtmParams);
                                [SST,f]=mtspectrumc((analogData(goodPos,STPos))',mtmParams);
                                plot(tfplotHandles(i,j),f,10*(log10(SST)-log10(SBL)),'color',plotColor); % in dB
                            else
                                disp('Choose same baseline and stimulus periods..');
                            end
                            if holdOnState
                                set(tfplotHandles(i,j),'Nextplot','add');
                            else
                                set(tfplotHandles(i,j),'Nextplot','replace');
                            end
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 1) % pcolor plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time
                            pcolor(tfplotHandles(i,j),t,f,log10(S1')); shading interp;
                            caxis([cmin cmax]);

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(log10(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(log10(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(i,j),t,f,dS1'); shading interp;
                            caxis([cmin cmax]);
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 2) % imagesc plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time
                            imagesc(t,f,log10(S1'),'Parent',tfplotHandles(i,j)); 
                            axis(tfplotHandles(i,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(i,j),'cLim',[cmin cmax]);

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(log10(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(log10(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum                       
                            imagesc(t,f,dS1','Parent',tfplotHandles(i,j)); 
                            axis(tfplotHandles(i,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(i,j),'cLim',[cmin cmax]);
                        end
                    end

                elseif (tfMethod == 2) % MP Algorithm
                    
                    disp(['number of atoms for reconstruction:' num2str(numAtomsMP)]);
                                       
                    if useBipolar
                        clear tag gaborInfo X
                        tag = 'temp/';
                        
                        disp('Running MP for bipolar case');
                        
                        % Import the data
                        X = analogData;
                        L = size(analogData,2);
                        signalRange = [1 L]; % full range
                        importData(X,mpFolderTemp,tag,signalRange,Fs);

                        % perform Gabor decomposition
                        Numb_points = L; % length of the signal
                        Max_iterations = 100; % number of iterations
                        
                        disp(['MP decomposition with max iterations =' num2str(Max_iterations)]);
                        
                        runGabor(folderName,tag,Numb_points, Max_iterations);
                        
                        gaborInfo = getGaborData(folderName,tag,1);
                        
                        gaborInfoGoodPos = gaborInfo{goodPos};
                        
                        wrap = [];
                        atomList = (1:numAtomsMP);
                        mpSpectrum = [];
                        disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                        for m=1:length(goodPos)
                            disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                            rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                            if m == 1
                                mpSpectrum = rEnergy;
                            else
                                mpSpectrum = mpSpectrum + rEnergy;
                            end
                        end
                        mpSpectrum = mpSpectrum/length(goodPos);
                        
                    else
                        clear tag gaborInfo
                        tag = channelString;
                        
                        load(fullfile(mpFolder,tag,'gaborInfo.mat'));
                        
                        gaborInfoGoodPos = gaborInfo(goodPos);
                        
                        L = size(analogData,2);
                        wrap = [];
                        atomList = (1:numAtomsMP);
                        mpSpectrum = [];
                        disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                        for m=1:length(goodPos)
                            disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                            rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                            if m == 1
                                mpSpectrum = rEnergy;
                            else
                                mpSpectrum = mpSpectrum + rEnergy;
                            end
                        end
                        mpSpectrum = mpSpectrum/length(goodPos);
                        
                    end
                    
                    t = timeVals;
                    f = 0:Fs/L:Fs/2; 
                    
                    % plot the MP Spectrum
                    if plotStyle == 3 % line plot
                        disp('No line plot for MP method');
                        
                    elseif plotStyle == 1 % pcolor plot
                        
                        if (spectrumType == 1) % raw
                            pcolor(tfplotHandles(i,j),t,f,log10(mpSpectrum)); 
                            shading interp;
                            caxis([cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(log10(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(log10(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(i,j),t,f,dmpSpectrum); 
                            shading interp;
                            caxis([cmin cmax]);
                        end
                        
                    elseif plotStyle == 2 % imagesc plot
                        
                        if (spectrumType == 1) % raw
                            imagesc(t,f,log10(mpSpectrum),'Parent',tfplotHandles(i,j)); 
                            axis(tfplotHandles(i,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(i,j),'cLim',[cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(log10(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(log10(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            imagesc(t,f,dmpSpectrum,'Parent',tfplotHandles(i,j)); 
                            axis(tfplotHandles(i,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(i,j),'cLim',[cmin cmax]);
                        end
                        
                    end
                    
                end


                % Display title
                if (i==1)
                    if (j==1)
                        title(tfplotHandles(i,j),['Azi: ' num2str(aValsUnique{1}(a))],'FontSize',titleFontSize);
                        % [Vinay] - because all the gabors are concentric,
                        % using the aValsUnique for the gabor 1 w.l.g.
                    else
                        title(tfplotHandles(i,j),num2str(aValsUnique{1}(a)),'FontSize',titleFontSize);
                    end
                end

                if (j==numCols)
                    if (i==1)
                     title(tfplotHandles(i,j),[{'Ele'} {num2str(eValsUnique{1}(e))}],'FontSize',titleFontSize,...
                         'Units','Normalized','Position',[1.25 0.5]);
                     % [Vinay] - because all the gabors are concentric,
                        % using the eValsUnique for the gabor 1 w.l.g.
                    else
                        title(tfplotHandles(i,j),num2str(eValsUnique{1}(e)),'FontSize',titleFontSize,...
                         'Units','Normalized','Position',[1.25 0.5]);
                    end
                end
            end
        end
    end
end

function tfplotLFPData1Parameter1Channel(tfplotHandles,channelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gaborNum,paramNum,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax,holdOnState)
        
        if ispc
            folderExtract = fullfile(folderName,'extractedData');
            folderSegment = fullfile(folderName,'segmentedData');
            mpFolder = fullfile(folderName,'mpAnalysis'); % for using stored MP data
            mpFolderTemp = 'C:\Documents\MP\data\'; % for online calculations
        else
            folderExtract = fullfile(folderName,'extractedData');
            folderSegment = fullfile(folderName,'segmentedData');
            mpFolder = fullfile(folderName,'mpAnalysis'); % for using stored MP data
            mpFolderTemp = '/home/vinay/Documents/MP/Data/'; % for online calculations
        end


        titleFontSize = 9;

%         timeForComputation = [40 100]/1000; % ms
%         freqForComputation = [40 60]; % Hz

        [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
            fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);
        [~,numCols] = size(tfplotHandles);

        % Get the data
        clear signal analogData
        load(fullfile(folderLFP,channelString));

        % Get bad trials
        badTrialFile = [folderSegment 'badTrials.mat'];
        if ~exist(badTrialFile,'file')
            disp('Bad trial file does not exist...');
            badTrials=[];
        else
            badTrials = loadBadTrials(badTrialFile);
            disp([num2str(length(badTrials)) ' bad trials']);
        end

        % [Vinay] - repeat the set of parameters for each gabor depending on the
        % number of columns to be drawn
        aList = repmat(a,1,numCols);
        eList = repmat(e,1,numCols);
        sList = repmat(s,1,numCols);
        fList = repmat(f,1,numCols);
        oList = repmat(o,1,numCols);
        cList = repmat(c,1,numCols);
        tList = repmat(t,1,numCols);
        rList = repmat(r,1,numCols);
        pList = repmat(p,1,numCols);

        % [Vinay] decide the parameter based on paramNum
        switch (paramNum-1)
            case 1
                aList(gaborNum,:) = 1:length(aValsUnique{gaborNum}); 
                titleParam = 'Azi: ';
                titleList = aValsUnique{gaborNum};
            case 2
                eList(gaborNum,:) = 1:length(eValsUnique{gaborNum}); 
                titleParam = 'Ele: ';
                titleList = eValsUnique{gaborNum};
            case 3
                sList(gaborNum,:) = 1:length(sValsUnique{gaborNum}); 
                titleParam = 'Sigma: ';
                titleList = sValsUnique{gaborNum};
            case 4
                fList(gaborNum,:) = 1:length(fValsUnique{gaborNum}); 
                titleParam = 'SF: ';
                titleList = fValsUnique{gaborNum};
            case 5
                oList(gaborNum,:) = 1:length(oValsUnique{gaborNum}); 
                titleParam = 'Ori: ';
                titleList = oValsUnique{gaborNum};
            case 6
                cList(gaborNum,:) = 1:length(cValsUnique{gaborNum}); 
                titleParam = 'Contr: ';
                titleList = cValsUnique{gaborNum};
            case 7
                tList(gaborNum,:) = 1:length(tValsUnique{gaborNum}); 
                titleParam = 'TF: ';
                titleList = tValsUnique{gaborNum};
            case 8
                rList(gaborNum,:) = 1:length(rValsUnique{gaborNum}); 
                titleParam = 'Rad: ';
                titleList = rValsUnique{gaborNum};
            case 9
                pList(gaborNum,:) = 1:length(pValsUnique{gaborNum}); 
                titleParam = 'SPhs: ';
                titleList = pValsUnique{gaborNum};
            otherwise
                disp('No particular gabor or parameter selected');
        end

        % [Vinay] - Get the lengths of indices in parameterCombinations
        for i = 1:3
            aLen(i) = length(aValsUnique{i});
            eLen(i) = length(eValsUnique{i});
            sLen(i) = length(sValsUnique{i});
            fLen(i) = length(fValsUnique{i});
            oLen(i) = length(oValsUnique{i});
            cLen(i) = length(cValsUnique{i});
            tLen(i) = length(tValsUnique{i});
            pLen(i) = length(pValsUnique{i}); % Vinay - for CRS
            rLen(i) = length(rValsUnique{i}); % Vinay - for CRS

            % If more than one value, then length is one greater for all the values
            % together
            if (aLen(i)> 1)           aLen(i)=aLen(i)+1;                    end
            if (eLen(i)> 1)           eLen(i)=eLen(i)+1;                    end
            if (sLen(i)> 1)           sLen(i)=sLen(i)+1;                    end
            if (fLen(i)> 1)           fLen(i)=fLen(i)+1;                    end
            if (oLen(i)> 1)           oLen(i)=oLen(i)+1;                    end
            if (cLen(i)> 1)           cLen(i)=cLen(i)+1;                    end
            if (tLen(i)> 1)           tLen(i)=tLen(i)+1;                    end
            if (pLen(i)> 1)           pLen(i)=pLen(i)+1;                    end % Vinay - for CRS
            if (rLen(i)> 1)           rLen(i)=rLen(i)+1;                    end % Vinay - for CRS
        end

        % Main loop
%         computationVals=zeros(1,numCols);
        for j=1:numCols
            clear goodPos
        %     goodPos = parameterCombinations{aList(j),eList(j),sList(j),fList(j),oList(j),cList(j),tList(j)};
        %     goodPos = setdiff(goodPos,badTrials);

            % [Vinay] - goodPos will be obtained by taking the goodPos for all
            % the three gabors and finding the common trials in them
            clear pos1 pos2 pos3
            switch protocolNumber
                       case 1 % Ring protocol
                           pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad is imp
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C

                       case 2 % Contrast Ring Protocol
                           pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2,j),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & cont
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C

                       case 3 % Dual Contrast Protocol
                           pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                           % [Vinay] - S is hidden in this case and only C & R
                           % define the stimuli. pos1 is therefore the full set
                           % here
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2,j),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & cont
                           % [Vinay] -  R params - cont & rad. Rest are matched to
                           % C's params
                           pos3 = parameterCombinations{aList(3,j),eList(3,j),sList(3,j),fList(3,j),oList(3,j),cList(3,j),tList(3,j),pList(3,j),rList(3,j),3}; % for C

                       case 4 % Dual Orientation Protocol
                           pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                           % [Vinay] - S is hidden in this case and only C & R
                           % define the stimuli. pos1 is therefore the full set
                           % here
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2,j),cLen(2),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & ori
                           % [Vinay] -  R params - ori & rad. Rest are matched to
                           % C's params
                           pos3 = parameterCombinations{aList(3,j),eList(3,j),sList(3,j),fList(3,j),oList(3,j),cList(3,j),tList(3,j),pList(3,j),rList(3,j),3}; % for C

                       case 5 % Dual Phase Protocol
                           pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                           % [Vinay] - S is hidden in this case and only C & R
                           % define the stimuli. pos1 is therefore the full set
                           % here
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2,j),rList(2,j),2}; % for R - only rad & p
                           % [Vinay] -  R params - p & rad. Rest are matched to
                           % C's params
                           pos3 = parameterCombinations{aList(3,j),eList(3,j),sList(3,j),fList(3,j),oList(3,j),cList(3,j),tList(3,j),pList(3,j),rList(3,j),3}; % for C

                       case 6 %Orientation Ring Protocol
                           pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2,j),cLen(2),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & ori
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C

                       case 7 %Phase Ring Protocol
                           pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2,j),rList(2,j),2}; % for R - only rad & p
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C

                       case 8 %Drifting Ring Protocol
                           pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tList(2,j),pLen(2),rList(2,j),2}; % for R - only rad & t
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rLen(3),3}; % for C

                       case 9 %Cross Orientation Protocol
                           disp('Left for future extension :)');
                       case 10 %Annulus Fixed Protocol
                           pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1}; % S params define C & S
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2,j),tLen(2),pLen(2),rList(2,j),2}; % for R - only rad & cont
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. However, the rad of C is
                           % still a relevant parameter and has to be taken care of
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3,j),3}; % for C - only rad

                otherwise
                    pos1 = parameterCombinations{aList(1,j),eList(1,j),sList(1,j),fList(1,j),oList(1,j),cList(1,j),tList(1,j),pList(1,j),rList(1,j),1};
                    % good positions for gabor 1 i.e. S
                    pos2 = parameterCombinations{aList(2,j),eList(2,j),sList(2,j),fList(2,j),oList(2,j),cList(2,j),tList(2,j),pList(2,j),rList(2,j),2}; % for R
                    pos3 = parameterCombinations{aList(3,j),eList(3,j),sList(3,j),fList(3,j),oList(3,j),cList(3,j),tList(3,j),pList(3,j),rList(3,j),3}; % for C
            end

            goodPos = intersect(pos1,pos2);
            goodPos = intersect(goodPos,pos3);
            goodPos = setdiff(goodPos,badTrials);

            if isempty(goodPos)
                disp('No entries for this combination..')
            else
                disp(['pos=' num2str(j) ',n=' num2str(length(goodPos))]);

                Fs = round(1/(timeVals(2)-timeVals(1)));
                BLRange = uint16((BLMax-BLMin)*Fs);
                STRange = uint16((STMax-STMin)*Fs);
                BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
                STPos = find(timeVals>=STMin,1)+ (1:STRange);

                xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
                xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);

        %         % Vinay - added next two 'if' conditions to take care of the case
        %         % when the lengths of BLPos(STPos) and xsBL(xsST) differ by 1
        %         if ((abs(length(BLRange)-length(xsBL)))==1)
        %             BLRange = ceil((BLMax-BLMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
        %         end
        %         
        %         if ((abs(length(STRange)-length(xsST)))==1)
        %             STRange = ceil((STMax-STMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
        %         end

%                 xsComputation = intersect(find(timeVals>=timeForComputation(1)),find(timeVals<timeForComputation(2)));
%                 freqComputation = intersect(find(xsST>=freqForComputation(1)),find(xsST<=freqForComputation(2)));

                % Vinay - added this notch data check
                if notchData
                    analogData = analogDataNotched;
                end

                if useBipolar
                    analogChannelString1 = channelString; % first electrode selected
                    clear analogData analogDataNotched
                    load(fullfile(folderLFP,analogChannelString1));

                    if notchData
                        analogData1 = analogDataNotched;
                    else
                        analogData1 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 1
                    end
                    % Vinay - store the analogData 
                    % for electrode 1 in a separate variable 
                    % otherwise loading the data for 2nd electrode will
                    % overwrite analogData

                    clear analogData analogDataNotched
                    load(fullfile(folderLFP,analogChannelString2));

                    if notchData
                        analogData2 = analogDataNotched;
                    else
                        analogData2 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 2
                    end
                    analogData = analogData1 - analogData2;
                end
                
                
                %--------------------------------------------------------------
                % Implement the method set by tfMethod varibale

                if (tfMethod == 1) % MTM
                    %----------------------------------------------------------
                    if (plotStyle == 3) % line plot (power vs f)
                        if (spectrumType == 1) % raw
                            [SBL,f]=mtspectrumc((analogData(goodPos,BLPos))',mtmParams); % baseline period
                            [SST,f]=mtspectrumc((analogData(goodPos,STPos))',mtmParams); % stimulus period
                            plot(tfplotHandles(j),f,log10(SST),'color',plotColor);
                            set(tfplotHandles(j),'Nextplot','add');
                            plot(tfplotHandles(j),f,log10(SBL),'color','g');
                            if holdOnState
                                set(tfplotHandles(j),'Nextplot','add');
                            else
                                set(tfplotHandles(j),'Nextplot','replace');
                            end

                        elseif (spectrumType == 2) % difference from baseline
                            if xsST == xsBL 
                                [SBL,f]=mtspectrumc((analogData(goodPos,BLPos))',mtmParams);
                                [SST,f]=mtspectrumc((analogData(goodPos,STPos))',mtmParams);
                                plot(tfplotHandles(j),f,10*(log10(SST)-log10(SBL)),'color',plotColor);
                            else
                                disp('Choose same baseline and stimulus periods..');
                            end
                            if holdOnState
                                set(tfplotHandles(j),'Nextplot','add');
                            else
                                set(tfplotHandles(j),'Nextplot','replace');
                            end
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 1) % pcolor plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time
                            pcolor(tfplotHandles(j),t,f,log10(S1')); shading interp;
                            caxis([cmin cmax]);

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(log10(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(log10(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(j),t,f,dS1'); shading interp;
                            caxis([cmin cmax]);
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 2) % imagesc plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time
                            imagesc(t,f,log10(S1'),'Parent',tfplotHandles(j)); 
                            axis(tfplotHandles(j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(j),'cLim',[cmin cmax]);

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(log10(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(log10(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum                       
                            imagesc(t,f,dS1','Parent',tfplotHandles(j)); 
                            axis(tfplotHandles(j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(j),'cLim',[cmin cmax]);
                        end
                    end
                    
                    %======================================================

                elseif (tfMethod == 2) % MP Algorithm
                    
                    disp(['number of atoms for reconstruction:' num2str(numAtomsMP)]);
                                       
                    if useBipolar
                        clear tag gaborInfo X
                        tag = 'temp/';
                        
                        disp('Running MP for the bipolar case...');
                        
                        % Import the data
                        X = analogData;
                        L = size(analogData,2);
                        signalRange = [1 L]; % full range
                        importData(X,mpFolderTemp,tag,signalRange,Fs);

                        % perform Gabor decomposition
                        Numb_points = L; % length of the signal
                        Max_iterations = 100; % number of iterations
                        
                        disp(['MP decomposition with max iterations =' num2str(Max_iterations)]);
                        
                        runGabor(folderName,tag,Numb_points, Max_iterations);
                        
                        gaborInfo = getGaborData(folderName,tag,1);
                        
                        gaborInfoGoodPos = gaborInfo{goodPos};
                        
                        wrap = [];
                        atomList = (1:numAtomsMP);
                        mpSpectrum = [];
                        disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                        for m=1:length(goodPos)
                            disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                            rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                            if m == 1
                                mpSpectrum = rEnergy;
                            else
                                mpSpectrum = mpSpectrum + rEnergy;
                            end
                        end
                        mpSpectrum = mpSpectrum/length(goodPos);
                        
                    else
                        clear tag gaborInfo
                        tag = channelString;
                        
                        load(fullfile(mpFolder,tag,'gaborInfo.mat'));
                        
                        gaborInfoGoodPos = gaborInfo(goodPos);
                        
                        L = size(analogData,2);
                        wrap = [];
                        atomList = (1:numAtomsMP);
                        mpSpectrum = [];
                        disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                        for m=1:length(goodPos)
                            disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                            rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                            if m == 1
                                mpSpectrum = rEnergy;
                            else
                                mpSpectrum = mpSpectrum + rEnergy;
                            end
                        end
                        mpSpectrum = mpSpectrum/length(goodPos);
                        
                    end
                    
                    t = timeVals;
                    f = 0:Fs/L:Fs/2; 
                    
                    % plot the MP Spectrum
                    if plotStyle == 3 % line plot
                        disp('No line plot for MP method');
                        
                    elseif plotStyle == 1 % pcolor plot
                        
                        if (spectrumType == 1) % raw
                            pcolor(tfplotHandles(j),t,f,log10(mpSpectrum)); 
                            shading interp;
                            caxis([cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(log10(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(log10(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(j),t,f,dmpSpectrum); 
                            shading interp;
                            caxis([cmin cmax]);
                        end
                        
                    elseif plotStyle == 2 % imagesc plot
                        
                        if (spectrumType == 1) % raw
                            imagesc(t,f,log10(mpSpectrum),'Parent',tfplotHandles(j)); 
                            axis(tfplotHandles(j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(j),'cLim',[cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(log10(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(log10(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            imagesc(t,f,dmpSpectrum,'Parent',tfplotHandles(j)); 
                            axis(tfplotHandles(j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(j),'cLim',[cmin cmax]);
                        end
                        
                    end
                    
                end
                
                % Display title
                if (j==1)
                    title(tfplotHandles(j),[titleParam num2str(titleList(j))],'FontSize',titleFontSize);
                else
                    title(tfplotHandles(j),num2str(titleList(j)),'FontSize',titleFontSize);
                end
                
                
            end
        end
        
end


