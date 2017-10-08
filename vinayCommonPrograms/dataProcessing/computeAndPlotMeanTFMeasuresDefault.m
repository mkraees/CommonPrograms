% plot mean TF across electrodes, sessions, conditions
% Vinay Shirhatti, 20 October 2016
% modified and developed from computeAndPlotMeanPSDMeasures
%
% INPUTS
% data : cell array of dimensions (numRows x numCols x numProtocols)
% numRows and numCols are number of different conditions (corresponding to
% certain stimulus parameters) varying along rows and columns respectively.
% numProtocols is the number of protocols to be analyzed together. Each
% cell is then itself a cell array of size (1 x numElecs) where numElecs is
% the number of electrodes to be analyzed
% eg: if there are 5 rows, 6 columns, 3 sessions and 12 electrodes then:
% data is (5 x 6 x 3) cell
% any data{r,c,i} is (1 x 12) cell array and each cell in this contains the
% trialwise signals for that electrode. For eg. if there are 16 trials for
% any particular condition for electrode number 3 in the list then:
% data{r,c,i}{3} : (16 x numTimeSamples) double array
%
% params : structure of variables/parameters that decide the computation
% drawplots : set this to 1 to draw the plots, otheriwse to 0 to only do
%             the calculations and return the TF data
% newfiguresflag: set this flag to 1 to plot everything. If the relevant 
%                 handles are passed then the plots are drawn in those
%                 subplots. Otherwise new figures are created and the plots
%                 are drawn in them.
% returnData : set this flag to 1 to return the relevant data, else tfdata
%              is empty
%
% OUTPUTS
% params : modified params
% psddata : the calculated TF data
% *************************************************************************
%
% This code is taken from computeAndPlotMeanTFMeasures. The computations
% here correspond to only the case of common baseline across conditions for
% every electrode and taking log at the level of every electrode before
% taking the mean across electrodes
%
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function [params,tfdata] = computeAndPlotMeanTFMeasuresDefault(data,params,drawplots,newfiguresflag,returnData)

if ~exist('drawplots','var'); drawplots=1; end
if ~exist('newfiguresflag','var'); newfiguresflag=0; end
if ~exist('returnData','var'); returnData=1; end

%--------------------------------------------------------------------------
% Read all sizes
numRows = size(data,1);
numCols = size(data,2);
numProtocols = size(data,3);
numElecs = size(data{1},2);

if ~isfield(params,'maxPlotsAlongX')
    params.maxPlotsAlongX = 6; % this variable stores the maximum number of plots to be drawn along the horizontal axis
end


% assign timeVals if not present
if ~isfield(params,'timeVals')
    params.timeVals = 1:size(data{1},1);
end

% convert the colors to a sequential list of colors if a conditionwise cell
% array has been passed
if strcmp(params.colorsListType,'matrix')
    tempColorsList = params.colorsList;
    params = rmfield(params,'colorsList');
    for r=1:numRows
        for c=1:numCols
            params.colorsList((r-1)*numCols+c,:) = tempColorsList{r,c};
        end
    end
end

if ~isfield(params,'mtmParams')
    params.mtmParams.tapers = [1 1]; % stft
    if strcmpi(params.gridType,'EEG')
        params.mtmParams.Fs = 1000;
    else
        params.mtmParams.Fs = 2000;
    end
    params.mtmParams.fpass = [params.fmin params.fmax];
    params.mtmParams.trialave = 0;
    params.mtmParams.err = 0;
    params.mtmParams.pad = -1; % this is essential to avoid padding the data while taking fft
else
    if ~isfield(params.mtmParams,'Fs')
        if strcmpi(params.gridType,'EEG')
            params.mtmParams.Fs = 1000;
        else
            params.mtmParams.Fs = 2000;
        end
    end
    if ~isfield(params.mtmParams,'fpass')
        params.mtmParams.fpass = [params.fmin params.fmax];
    end
    if ~isfield(params.mtmParams,'trialave')
        params.mtmParams.trialave = 0;
    end
    if ~isfield(params.mtmParams,'err')
        params.mtmParams.err = [1 0.05];
    end
    if ~isfield(params.mtmParams,'pad')
        params.mtmParams.pad = -1;
    end
end

if ~isfield(params,'movingWin')
    params.movingWin = [0.25 0.025]; % default
end

%--------------------------------------------------------------------------
% Compute and plot TFs
%--------------------------------------------------------------------------
% Initialize variables
TFEveryTrial = cell(numRows,numCols,numProtocols);
TFAcrossTrials = cell(numRows,numCols,numProtocols);
serrTFAcrossTrials = cell(numRows,numCols,numProtocols);

TFAcrossTrialsBL = cell(numRows,numCols,numProtocols);

meanTFAcrossTrialsLogElectrodes = cell(numRows,numCols,numProtocols);

% Compute TF and related measures
for i = 1:numProtocols
    numElecs = size(data{1,1,i},2);
    for r = 1:numRows
        for c = 1:numCols
            disp([' TF >> processing: protocol ' num2str(i) '/' num2str(numProtocols) ', row ' num2str(c) '/' num2str(numCols) ', col ' num2str(r) '/' num2str(numRows)]);
            
            clear dataBL dataST
            
            [TFEveryTrial{r,c,i},t,f,serrTFAcrossTrials{r,c,i}] = cellfun(@(y) mtspecgramc(y,params.movingWin,params.mtmParams),cellfun(@(x) x',data{r,c,i},'UniformOutput',0),'UniformOutput',0);
            % cellfun(@(x) x',data{r,c,i},'UniformOutput',0) takes the
            % conditionwise signals for every electrode, every trial and
            % transposes them so that they are of the form: samples x
            % trials, as required by mtspectrumc. Thus the result of this
            % operation is a cell array (of size 1 x numElecs) consisting
            % of signals from all electrodes with each cell corresponding
            % to all trials from that electrode. mtspectrumc is then run on
            % each electrode's trials. Looks great! cellfun's lot of fun!
            % :D
            
            TFAcrossTrials{r,c,i} = cellfun(@(x) squeeze(nanmean(x,3)),TFEveryTrial{r,c,i},'UniformOutput',0);
            % mean across trials for every electrode => mean TF for each
            % electrode for this condition
            
            params.TFtimeVals = params.timeVals(1) + t{1}; 
            params.freqVals = f{1};
            numTimePoints = length(params.TFtimeVals);
            numFreqPoints = length(params.freqVals);
            
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common')
                % weighted mean: the number of repeats per electrode may 
                % vary. In such a case one may take a weighted mean 
                % This approach is required only when you have unequal
                % repeats per electrode i.e. when the badTrials selection
                % strategy is not the 'common' badTrials strategy
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:))); 
                meanTFAcrossTrialsLogElectrodes{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) conv2Log(x)*y./sum(weights), TFAcrossTrials{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
                % cell2mat(cellfun(@(x,y) x*y./sum(weights),
                % TFAcrossTrials{r,c,i}, num2cell(weights'),'UniformOutput',0);
                % This multiplies the mean spectrum for each electrode with
                % the correspoding number of repeats used to generate it
                % and divides by the total number of repeats summed across
                % all electrodes => weighted contribution of each electrode
                % This is converted to a matrix (i.e. cell2mat) and
                % reshaped so that the dimensions are time x freq x
                % numElecs. We then sum across all the electrodes to get
                % the mean TF spectrum across trials and electrodes for
                % this condition and session i.e. squeeze(nansum(x,3))

            else
                meanTFAcrossTrialsLogElectrodes{r,c,i} = squeeze(nanmean(reshape(cell2mat([TFAcrossTrials{r,c,i}]),numTimePoints,numFreqPoints,numElecs),3));
                % Take spectra from all electrodes together: [x{r,c,i}],
                % convert it to a matrix (cell2mat(x)), reshape it to
                % dimensions time x freq x numElecs, take mean across the
                % electrodes and squeeze i.e. squeeze(nanmean(x,3))
            end
            
            %--------------------------------------------------------------
            % Baseline calculation
            % get the baseline and stimulus indices
            blPos = params.TFtimeVals>=params.blmin & params.TFtimeVals<=params.blmax;
            stPos = params.TFtimeVals>=params.stmin & params.TFtimeVals<=params.stmax;
            
            % mean across the baseline time indices to give the mean power at every frequency in the baseline period 
            TFAcrossTrialsBL{r,c,i} = cellfun(@(x) mean(x(blPos,:),1), TFAcrossTrials{r,c,i}, 'UniformOutput',0);
        end
    end
end

%--------------------------------------------------------------------------
% Baseline calculation
%--------------------------------------------------------------------------

% Compute common baseline across conditions per electrode
commonBaselineTFAcrossTrials = cell(1,numProtocols);
TFAcrossTrialsBLTranspose = cellfun(@(x) x',TFAcrossTrialsBL,'UniformOutput',0); % to align the electrodes along the rows
for i = 1:numProtocols % find the common baseline for each session separately
    numElecs = size(data{1,1,i},2);
    commonBaselineTFAcrossTrials{i} = cell(numElecs,1);
    TFAcrossTrialsAllConditions = [TFAcrossTrialsBLTranspose{:,:,i}];
    for en = 1:numElecs
        commonBaselineTFAcrossTrials{i}{en,:} = mean(cell2mat(TFAcrossTrialsAllConditions(en,:)'),1);
    end
    % take all conditions together for a protocol and find the mean
    % baseline TF across all conditions for each electrode separately.
    % This is done by taking the mean baseline across conditions
end


%--------------------------------------------------------------------------
% Calculate TF measures: difference TF spectrum, band power time course,
% peak frequency time time course
%--------------------------------------------------------------------------
dTFAcrossTrialsCommonBL = cell(numRows,numCols,numProtocols);
dTFAcrossTrialsLogElectrodesCommonBL = cell(numRows,numCols,numProtocols);

dBandPowerAcrossTrialsCommonBL = cell(numRows,numCols,numProtocols);
dBandPowerAcrossTrialsLogElectrodesCommonBL = cell(numRows,numCols,numProtocols);

peakFreqdBandPowerAcrossTrials = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossTrialsCommonBL = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);

for r = 1:numRows
    for c = 1:numCols
        for i = 1:numProtocols
            numElecs = size(data{1,1,i},2);
            %__Difference TFs(log) wrt common baseline across conditions__
            clear commonBLTFAcrossTrials commonBLTFAcrossLogTrials
            % repeat the mean baseline power per frequency across all time
            % points and subtract it from the spectrum
            commonBLAcrossTrials = cellfun(@(x) repmat(x,numTimePoints,1), commonBaselineTFAcrossTrials{i}, 'UniformOutput',0);
            dTFAcrossTrialsCommonBL{r,c,i} = cellfun(@(x,y) conv2Log(x./y),TFAcrossTrials{r,c,i},commonBLAcrossTrials','UniformOutput',0);
            
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common') % averaging across electrodes with possibly unequal repeats
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:))); % number of repeats per electrode for a session
                dTFAcrossTrialsLogElectrodesCommonBL{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) x*y./sum(weights), dTFAcrossTrialsCommonBL{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
                % Multiply the mean diff spectrum for each electrode with
                % the correspoding number of repeats used to generate it
                % and divide by the total number of repeats summed across
                % all electrodes => weighted contribution of each electrode
                % This is converted to a matrix (i.e. cell2mat) and
                % reshaped so that the dimensions are time x freq x
                % numElecs. We then sum across all the electrodes to get
                % the mean TF spectrum across trials and electrodes for
                % this condition and session i.e. squeeze(nansum(x,3))
            else
                dTFAcrossTrialsLogElectrodesCommonBL{r,c,i} = squeeze(nanmean(reshape(cell2mat([dTFAcrossTrialsCommonBL{r,c,i}]),numTimePoints,numFreqPoints,numElecs),3));
            end

            %____________________TF Measures_______________________
            fPos = params.freqVals>=params.fBandLow & params.freqVals<=params.fBandHigh; % frequency indices for the band of interest
            if isfield(params,'extraFreqBands')
                for l=1:length(params.extraFreqBands)
                    fPosExtra{l} = params.freqVals>=params.extraFreqBands{l}(1) & params.freqVals<=params.extraFreqBands{l}(2);
                end
            end
            
            % Power in the selected frequency band
            % Time course of the difference power in the selected band
            dBandPowerAcrossTrialsCommonBL{r,c,i} = cellfun(@(x,y) conv2Log(sum(x(:,fPos),2)./sum(y(:,fPos),2)),TFAcrossTrials{r,c,i},commonBLAcrossTrials','UniformOutput',0);
            
            % Averaged across electrodes
            dBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,i} = nanmean(cell2mat(dBandPowerAcrossTrialsCommonBL{r,c,i}),2);
            
            % Peak frequency in the selected frequency band
            bandFreqVals = params.freqVals(fPos);
            
            clear peakPowerCommonBL peakPowerLogCommonBL

            peakPowerCommonBL = cellfun(@(x) max(x(:,fPos),[],2), dTFAcrossTrialsCommonBL{r,c,i}, 'uniformoutput',0);
            for tx=1:numTimePoints
                peakFreqdBandPowerAcrossTrialsCommonBL{r,c,i}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPos)==y(tx)), dTFAcrossTrialsCommonBL{r,c,i}, peakPowerCommonBL, 'uniformoutput',0));
            end
            
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common')
                clear weights
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:))); % number of repeats per electrode for a session
                peakFreqdBandPowerAcrossTrialsElectrodesCommonBL{r,c,i} = peakFreqdBandPowerAcrossTrialsCommonBL{r,c,i}*weights./sum(weights);
            else
                peakFreqdBandPowerAcrossTrialsElectrodesCommonBL{r,c,i} = mean(peakFreqdBandPowerAcrossTrialsCommonBL{r,c,i},2);
            end
        end
    end
end

%==========================================================================
% Computations for the extra frequency bands, if specified
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isfield(params,'extraFreqBands')
    
    for l=1:length(params.extraFreqBands)
        fPosExtra{l} = params.freqVals>=params.extraFreqBands{l}(1) & params.freqVals<=params.extraFreqBands{l}(2);
    end
    
    numBands = length(params.extraFreqBands);

    dBandPowerAcrossTrialsCommonBLEB = cell(numRows,numCols,numProtocols,numBands);
    
    dBandPowerAcrossTrialsLogElectrodesCommonBLEB = cell(numRows,numCols,numProtocols,numBands);

    peakFreqdBandPowerAcrossTrialsCommonBLEB = cell(numRows,numCols,numProtocols,numBands);

    peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB = cell(numRows,numCols,numProtocols,numBands);

    for r = 1:numRows
        for c = 1:numCols
            for i = 1:numProtocols
                numElecs = size(data{1,1,i},2);
                commonBLAcrossTrials = cellfun(@(x) repmat(x,numTimePoints,1), commonBaselineTFAcrossTrials{i}, 'UniformOutput',0);
                for l=1:numBands
                    
                    % Power in the selected frequency band
                    % Time course of the difference power in the selected band
                    
                    dBandPowerAcrossTrialsCommonBLEB{r,c,i,l} = cellfun(@(x,y) conv2Log(sum(x(:,fPosExtra{l}),2)./sum(y(:,fPosExtra{l}),2)),TFAcrossTrials{r,c,i},commonBLAcrossTrials','UniformOutput',0);
                    
                    % Averaged across electrodes
                    dBandPowerAcrossTrialsLogElectrodesCommonBLEB{r,c,i,l} = nanmean(cell2mat(dBandPowerAcrossTrialsCommonBLEB{r,c,i,l}),2);
                    
                    % Peak frequency in the selected frequency band
                    bandFreqVals = params.freqVals(fPosExtra{l});

                    clear peakPowerCommonBL
                    peakPowerCommonBL = cellfun(@(x) max(x(:,fPosExtra{l}),[],2), dTFAcrossTrialsCommonBL{r,c,i}, 'uniformoutput',0);
                    for tx=1:numTimePoints
                        peakFreqdBandPowerAcrossTrialsCommonBLEB{r,c,i,l}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPosExtra{l})==y(tx)), dTFAcrossTrialsCommonBL{r,c,i}, peakPowerCommonBL, 'uniformoutput',0));
                    end

                    if params.weightedmean && ~strcmpi(params.badtrialsoption,'common')
                        clear weights
                        weights = cell2mat(squeeze(params.numRepeats(r,c,i,:))); % number of repeats per electrode for a session
                        peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB{r,c,i,l} = peakFreqdBandPowerAcrossTrialsCommonBL{r,c,i}*weights./sum(weights);
                    else
                        peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB{r,c,i,l} = mean(peakFreqdBandPowerAcrossTrials{r,c,i},2);
                    end   
                end
            end
        end
    end
end


if drawplots
%==========================================================================
%--------------------------------------------------------------------------
% Draw the plots
%--------------------------------------------------------------------------
% Plot TFs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Difference TFs: change of power from baseline to stimulus period
% With common baseline across conditions per electrode per session

% Check if plot handles have been passed, else create a new figure and plot
% handles if newfiguresflag is set
if ~isfield(params,'hdiffTFCommonBLFig') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hdiffTFCommonBLFig = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hdiffTFCommonBLFig')
    dTFAcrossTrialsLogElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints,numFreqPoints);
    meanTFAcrossTrialsLogElectrodesSessions = zeros(numRows,numCols,numTimePoints,numFreqPoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:))),2)';
                    dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dTFAcrossTrialsLogElectrodesCommonBL(r,c,:),1,[]), num2cell(weights),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the total dTF spectrum
                    meanTFAcrossTrialsLogElectrodesSessions(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(meanTFAcrossTrialsLogElectrodes(r,c,:),1,[]), num2cell(weights),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                else
                    dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:) = squeeze(nanmean(cell2mat(dTFAcrossTrialsLogElectrodesCommonBL(r,c,:)),3));
                    meanTFAcrossTrialsLogElectrodesSessions(r,c,:,:) = squeeze(nanmean(cell2mat(meanTFAcrossTrialsLogElectrodes(r,c,:)),3));
                end
            else
                dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:) = dTFAcrossTrialsLogElectrodesCommonBL{r,c,:};
                meanTFAcrossTrialsLogElectrodesSessions(r,c,:,:) = meanTFAcrossTrialsLogElectrodes{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hdiffTFCommonBLFig(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hdiffTFCommonBLFig(x,y));
            end
%             pcolor(params.TFtimeVals,params.freqVals,10*squeeze(dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:))');
%             pcolorAsImage(params.hdiffTFCommonBLFig(x,y),params.TFtimeVals,params.freqVals,10*squeeze(dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:))');
%             setAxesProperties(params.hdiffTFCommonBLFig(x,y),params,x,y,rows);
%             showTitleAndN(params.hdiffTFCommonBLFig(x,y),params,r,c);
%             shading interp; colormap(jet);
            tfXLims = [params.tmin params.tmax];
            tfYLims = [params.fmin params.fmax];
            cLims = [params.cmin params.cmax];
            CData = 10*squeeze(dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:))';
            spHandle = subplot(params.hdiffTFCommonBLFig(x,y)); 
            [imData] = plotPColorAsImage(spHandle,params.TFtimeVals,params.freqVals,CData,tfXLims,tfYLims,cLims);

            timeAxisTF = tfXLims;
            freqAxisTF = [params.fmin params.fmax];
            imagesc(timeAxisTF,fliplr(freqAxisTF),imData);
            set(params.hdiffTFCommonBLFig(x,y),'ydir','normal');
            
            setAxesProperties(params.hdiffTFCommonBLFig(x,y),params,x,y,rows);
            showTitleAndN(params.hdiffTFCommonBLFig(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end
    set(params.hdiffTFCommonBLFig,'clim',[params.cmin params.cmax]);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean diff TF across trials, electrodes and sessions (common baseline)','FitBoxToText','on');
    end
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hdiffTFCommonBLFig(r,c));
                delete(axisHandle);
            end
        end
    end
end

%==========================================================================
%--------------------------------------------------------------------------
% Draw the time course of band power and peak frequency
%--------------------------------------------------------------------------
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% with common baseline

if ~isfield(params,'hBandPowerWithTimeCommonBL') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hBandPowerWithTimeCommonBL = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hBandPowerWithTimeCommonBL')
    dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:))),2)';
                    dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dBandPowerAcrossTrialsLogElectrodesCommonBL(r,c,:),1,[]), num2cell(weights),'uniformoutput',0)),2));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the total dBandPower
                else
                    dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(cell2mat(dBandPowerAcrossTrialsLogElectrodesCommonBL(r,c,:)),3));
                end
            else
                dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = dBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hBandPowerWithTimeCommonBL(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hBandPowerWithTimeCommonBL(x,y));
            end
            plot(params.TFtimeVals,10*squeeze(dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
            hold on;
            line([params.timeVals(1) params.timeVals(end)],[0 0],'color',[0.5 0.4 0.7]);
            setAxesProperties(params.hBandPowerWithTimeCommonBL(x,y),params,x,y,rows);
            showTitleAndN(params.hBandPowerWithTimeCommonBL(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end

    set(params.hBandPowerWithTimeCommonBL,'ylim',ylimits);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes{1})],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean TF across trials, electrodes and sessions','FitBoxToText','on');
    end
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hBandPowerWithTimeCommonBL(r,c));
                delete(axisHandle);
            end
        end
    end
end

%==========================================================================
% peak frequency time course
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% with common baseline
if ~isfield(params,'hPeakFreqWithTimeCommonBL') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hPeakFreqWithTimeCommonBL = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hPeakFreqWithTimeCommonBL')
    peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:))),2)';
                    peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(peakFreqdBandPowerAcrossTrialsElectrodesCommonBL(r,c,:),1,[]), num2cell(weights),'uniformoutput',0)),2));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the overall peak freq
                else
                    peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(cell2mat(peakFreqdBandPowerAcrossTrialsElectrodesCommonBL(r,c,:)),3));
                end
            else
                peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = peakFreqdBandPowerAcrossTrialsElectrodesCommonBL{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hPeakFreqWithTimeCommonBL(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hPeakFreqWithTimeCommonBL(x,y));
            end
            plot(params.TFtimeVals,squeeze(peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
            hold on;
            line([params.timeVals(1) params.timeVals(end)],[0 0],'color',[0.5 0.4 0.7]);
            setAxesProperties(params.hPeakFreqWithTimeCommonBL(x,y),params,x,y,rows);
            showTitleAndN(params.hPeakFreqWithTimeCommonBL(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end

    set(params.hPeakFreqWithTimeCommonBL,'ylim',[params.fBandLow params.fBandHigh]);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean TF across trials, electrodes and sessions','FitBoxToText','on');
    end
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hPeakFreqWithTimeCommonBL(r,c));
                delete(axisHandle);
            end
        end
    end
end

%--------------------------------------------------------------------------
end
%==========================================================================
%==========================================================================
% return the specified tf data
%--------------------------------------------------------------------------
if returnData
    % default - take log at the level of every electrode
    if params.returnAllTF
        tfdata.TF = TFAcrossTrials;
        tfdata.TFMean = meanTFAcrossTrialsLogElectrodes;
        tfdata.diffTF = dTFAcrossTrialsCommonBL;
        tfdata.commonBaselineTFAcrossTrials = commonBaselineTFAcrossTrials;
    end
    tfdata.diffTFMean = dTFAcrossTrialsLogElectrodesCommonBL;
    tfdata.dBandPower = dBandPowerAcrossTrialsCommonBL;
    tfdata.dBandPowerMean = dBandPowerAcrossTrialsLogElectrodesCommonBL;
    tfdata.peakFreq = peakFreqdBandPowerAcrossTrialsCommonBL;
    tfdata.peakFreqMean = peakFreqdBandPowerAcrossTrialsElectrodesCommonBL;
    if drawplots && newfiguresflag && params.returnAllTF
        tfdata.TFMeanSessions = meanTFAcrossTrialsLogElectrodesSessions;
        tfdata.diffTFMeanSessions = dTFAcrossTrialsLogElectrodesSessionsCommonBL;
        tfdata.dBandPowerMeanSessions = dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL;
        tfdata.peakFreqMeanSessions = peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL;
    end
    if isfield(params,'extraFreqBands')
        tfdata.dBandPowerEB = dBandPowerAcrossTrialsCommonBLEB;
        tfdata.dBandPowerMeanEB = dBandPowerAcrossTrialsLogElectrodesCommonBLEB;
        tfdata.peakFreqEB = peakFreqdBandPowerAcrossTrialsCommonBLEB;
        tfdata.peakFreqMeanEB = peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB;
    end
    tfdata.TFtimeVals = params.TFtimeVals;
    tfdata.TFfreqVals = params.freqVals;
    tfdata.params = params;
else
    tfdata = [];
end


end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%==========================================================================
% Additional Functions
%==========================================================================

function setAxesProperties(hAxes,params,row,col,rows)

if ~isfield(params,'setAxesProperties')
    params.setAxesProperties = 1;
end

if params.setAxesProperties
    set(hAxes,'xlim',[params.tmin+params.movingWin(1)/2 params.tmax-params.movingWin(1)/2]);
    set(hAxes,'tickLength',[0.03 0.01]);
    set(hAxes,'fontsize',20,'fontweight','bold');
    set(hAxes,'xtick',[0 0.8],'tickdir','out');
    set(hAxes,'tickdir','out');
    set(hAxes,'box','off');

    if col==1 && row==rows
        ylabel(hAxes,'Frequency (Hz)');
        xlabel(hAxes,'Time (s)');
    else
        set(hAxes,'yticklabel',[]);
        set(hAxes,'xticklabel',[]);
    end
elseif params.cleanAxes
    set(hAxes,'yticklabel',[]);
    set(hAxes,'xticklabel',[]);
else
    setAxesBasicProperties(hAxes,params);
end
end

function showTitleAndN(hAxes,params,row,col,showLegend,x,y)

if ~exist('showLegend','var'); showLegend = 0; end
if ~exist('x','var'); x = 0.1; end
if ~exist('y','var'); y = 0.94; end

if params.showtitlex && params.showtitley
    if showLegend
        legend(hAxes,[params.titleStringX{:} ',' params.titleStringY{:}]);
    else
        text(x, y, [params.titleStringX{row,col} ',' params.titleStringY{row,col}],'fontsize',6,'unit','normalized');
    end
elseif params.showtitlex && ~params.showtitley
    if showLegend
        legend(hAxes,params.titleStringX{:},'fontsize',6);
    else
        text(x, y, params.titleStringX{row,col},'unit','normalized','fontsize',6);
    end
elseif ~params.showtitlex && params.showtitley
    if showLegend
        legend(hAxes,params.titleStringY{:});
    else
        text(x, y, params.titleStringY{row,col},'unit','normalized','fontsize',6);
    end
end

if params.showN
    meanN = mean(mean(squeeze(cell2mat(params.numRepeats(row,col,:,:)))'));
    text(x, y-0.1, ['n=' num2str(meanN)],'unit','normalized','fontsize',6);
end

end

function setAxesBasicProperties(hAxes,params)
set(hAxes,'xlim',[params.tmin+params.movingWin(1)/2 params.tmax-params.movingWin(1)/2]);
set(hAxes,'tickLength',[0.03 0.01]);
set(hAxes,'fontsize',20,'fontweight','bold');
set(hAxes,'tickdir','out');
set(hAxes,'box','off');  
end

function setAxesPropertiesFigure(hFig)
allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'xlim',[-0.5 1]);
set(allAxesInFigure,'tickLength',[0.03 0.01]);
set(allAxesInFigure,'fontsize',20,'fontweight','bold');
set(allAxesInFigure,'xtick',[0 0.8],'tickdir','out');
set(allAxesInFigure,'box','off');
end