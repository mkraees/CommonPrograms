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

function [params,tfdata] = computeAndPlotMeanTFMeasures(data,params,drawplots,newfiguresflag,returnData)

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
TFAcrossLogTrials = cell(numRows,numCols,numProtocols);

meanTFAcrossTrialsElectrodes = cell(numRows,numCols,numProtocols);
meanTFAcrossTrialsLogElectrodes = cell(numRows,numCols,numProtocols);
meanTFAcrossLogTrialsElectrodes = cell(numRows,numCols,numProtocols);

TFAcrossTrialsBL = cell(numElecs,numProtocols);
TFAcrossLogTrialsBL = cell(numElecs,numProtocols);

% Compute TF and related measures
for i = 1:numProtocols
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
            
            % Similarly mean across logarithm of every trial
            TFAcrossLogTrials{r,c,i} = cellfun(@(x) squeeze(nanmean(conv2Log(x),3)),TFEveryTrial{r,c,i},'UniformOutput',0);
            
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
                meanTFAcrossTrialsElectrodes{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) conv2Log(x*y./sum(weights)), TFAcrossTrials{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
                meanTFAcrossTrialsLogElectrodes{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) conv2Log(x)*y./sum(weights), TFAcrossTrials{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
                meanTFAcrossLogTrialsElectrodes{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) x*y./sum(weights), TFAcrossLogTrials{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
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
                meanTFAcrossTrialsElectrodes{r,c,i} = squeeze(conv2Log(nanmean(reshape(10.^cell2mat([TFAcrossTrials{r,c,i}]),numTimePoints,numFreqPoints,numElecs),3)));
                meanTFAcrossTrialsLogElectrodes{r,c,i} = squeeze(nanmean(reshape(cell2mat([TFAcrossTrials{r,c,i}]),numTimePoints,numFreqPoints,numElecs),3));
                meanTFAcrossLogTrialsElectrodes{r,c,i} = squeeze(nanmean(reshape(cell2mat([TFAcrossLogTrials{r,c,i}]),numTimePoints,numFreqPoints,numElecs),3));
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
            TFAcrossLogTrialsBL{r,c,i} = cellfun(@(x) mean(x(blPos,:),1), TFAcrossLogTrials{r,c,i}, 'UniformOutput',0);
            
        end
    end
end

%--------------------------------------------------------------------------
% Baseline calculation
%--------------------------------------------------------------------------

% Compute common baseline across conditions per electrode
commonBaselineTFAcrossTrials = cell(numElecs,numProtocols);
commonBaselineTFAcrossLogTrials = cell(numElecs,numProtocols);

TFAcrossTrialsBLTranspose = cellfun(@(x) x',TFAcrossTrialsBL,'UniformOutput',0); % to align the electrodes along the rows
TFAcrossLogTrialsBLTranspose = cellfun(@(x) x',TFAcrossLogTrialsBL,'UniformOutput',0); % to align the electrodes along the rows

for i = 1:numProtocols % find th common baseline for each session separately
    TFAcrossTrialsAllConditions = [TFAcrossTrialsBLTranspose{:,:,i}];
    TFAcrossLogTrialsAllConditions = [TFAcrossLogTrialsBLTranspose{:,:,i}];
    for en = 1:numElecs
        commonBaselineTFAcrossTrials{en,i} = mean(cell2mat(TFAcrossTrialsAllConditions(en,:)'),1);
        commonBaselineTFAcrossLogTrials{en,i} = mean(cell2mat(TFAcrossLogTrialsAllConditions(en,:)'),1);
    end
    % take all conditions together for a protocol and find the mean
    % baseline TF across all conditions for each electrode separately.
    % This is done by taking the mean baseline across conditions
end


%--------------------------------------------------------------------------
% Calculate TF measures: difference TF spectrum, band power time course,
% peak frequency time time course
%--------------------------------------------------------------------------

dTFAcrossTrials = cell(numRows,numCols,numProtocols);
dTFAcrossLogTrials = cell(numRows,numCols,numProtocols);

dTFAcrossTrialsElectrodes = cell(numRows,numCols,numProtocols);
dTFAcrossTrialsLogElectrodes = cell(numRows,numCols,numProtocols);
dTFAcrossLogTrialsElectrodes = cell(numRows,numCols,numProtocols);

dTFAcrossTrialsCommonBL = cell(numRows,numCols,numProtocols);
dTFAcrossLogTrialsCommonBL = cell(numRows,numCols,numProtocols);

dTFAcrossTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);
dTFAcrossTrialsLogElectrodesCommonBL = cell(numRows,numCols,numProtocols);
dTFAcrossLogTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);

dBandPowerAcrossTrials = cell(numRows,numCols,numProtocols);
dBandPowerAcrossLogTrials = cell(numRows,numCols,numProtocols);

dBandPowerAcrossTrialsCommonBL = cell(numRows,numCols,numProtocols);
dBandPowerAcrossLogTrialsCommonBL = cell(numRows,numCols,numProtocols);

dBandPowerAcrossTrialsElectrodes = cell(numRows,numCols,numProtocols);
dBandPowerAcrossTrialsLogElectrodes = cell(numRows,numCols,numProtocols);
dBandPowerAcrossLogTrialsElectrodes = cell(numRows,numCols,numProtocols);

dBandPowerAcrossTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);
dBandPowerAcrossTrialsLogElectrodesCommonBL = cell(numRows,numCols,numProtocols);
dBandPowerAcrossLogTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);

peakFreqdBandPowerAcrossTrials = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossLogTrials= cell(numRows,numCols,numProtocols);

peakFreqdBandPowerAcrossTrialsCommonBL = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossLogTrialsCommonBL = cell(numRows,numCols,numProtocols);

peakFreqdBandPowerAcrossTrialsElectrodes = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossLogTrialsElectrodes = cell(numRows,numCols,numProtocols);

peakFreqdBandPowerAcrossTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);

for r = 1:numRows
    for c = 1:numCols
        for i = 1:numProtocols
            %______________Difference TFs (logarithmic)___________________
            % per electrode averaged across its trials
            clear baselineTFAcrossTrials baselineTFAcrossLogTrials
            % repeat the mean baseline power per frequency across all time
            % points and subtract it from the spectrum
            baselineTFAcrossTrials = cellfun(@(x) repmat(x,numTimePoints,1), TFAcrossTrialsBL{r,c,i}, 'UniformOutput',0);
            baselineTFAcrossLogTrials = cellfun(@(x) repmat(x,numTimePoints,1), TFAcrossTrialsBL{r,c,i}, 'UniformOutput',0);
            dTFAcrossTrials{r,c,i} = cellfun(@(x,y) conv2Log(x./y),TFAcrossTrials{r,c,i},baselineTFAcrossTrials,'UniformOutput',0);
            dTFAcrossLogTrials{r,c,i} = cellfun(@(x,y) (x-y),TFAcrossLogTrials{r,c,i},baselineTFAcrossLogTrials,'UniformOutput',0);
            
            
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common') % averaging across electrodes with possibly unequal repeats
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:))); % number of repeats per electrode for a session
                dTFAcrossTrialsElectrodes{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) conv2Log((10.^x)*y./sum(weights)), dTFAcrossTrials{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
                dTFAcrossTrialsLogElectrodes{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) x*y./sum(weights), dTFAcrossTrials{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
                dTFAcrossLogTrialsElectrodes{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) x*y./sum(weights), dTFAcrossLogTrials{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
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
                % calculate mean across electrodes
                dTFAcrossTrialsElectrodes{r,c,i} = squeeze(conv2Log(mean(reshape(10.^cell2mat(dTFAcrossTrials{r,c,i}),numTimePoints,numFreqPoints,numElecs),3)));
                dTFAcrossTrialsLogElectrodes{r,c,i} = squeeze(mean(reshape(cell2mat(dTFAcrossTrials{r,c,i}),numTimePoints,numFreqPoints,numElecs),3));
                dTFAcrossLogTrialsElectrodes{r,c,i} = squeeze(mean(reshape(cell2mat(dTFAcrossLogTrials{r,c,i}),numTimePoints,numFreqPoints,numElecs),3));
            end
            
            %__Difference TFs(log) wrt common baseline across conditions__
            clear commonBLTFAcrossTrials commonBLTFAcrossLogTrials
            % repeat the mean baseline power per frequency across all time
            % points and subtract it from the spectrum
            commonBLAcrossTrials = cellfun(@(x) repmat(x,numTimePoints,1), commonBaselineTFAcrossTrials(:,i), 'UniformOutput',0);
            commonBLAcrossLogTrials = cellfun(@(x) repmat(x,numTimePoints,1), commonBaselineTFAcrossLogTrials(:,i), 'UniformOutput',0);
            dTFAcrossTrialsCommonBL{r,c,i} = cellfun(@(x,y) conv2Log(x./y),TFAcrossTrials{r,c,i},commonBLAcrossTrials','UniformOutput',0);
            dTFAcrossLogTrialsCommonBL{r,c,i} = cellfun(@(x,y) (x-y),TFAcrossLogTrials{r,c,i},commonBLAcrossLogTrials','UniformOutput',0);
            
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common') % averaging across electrodes with possibly unequal repeats
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:))); % number of repeats per electrode for a session
                dTFAcrossTrialsElectrodesCommonBL{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) conv2Log((10.^x)*y./sum(weights)), dTFAcrossTrialsCommonBL{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
                dTFAcrossTrialsLogElectrodesCommonBL{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) x*y./sum(weights), dTFAcrossTrialsCommonBL{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
                dTFAcrossLogTrialsElectrodesCommonBL{r,c,i} = squeeze(nansum(reshape(cell2mat(cellfun(@(x,y) x*y./sum(weights), dTFAcrossLogTrialsCommonBL{r,c,i}, num2cell(weights'),'UniformOutput',0)),numTimePoints,numFreqPoints,numElecs),3));
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
                dTFAcrossTrialsElectrodesCommonBL{r,c,i} = squeeze(conv2Log(nanmean(reshape(10.^cell2mat([dTFAcrossTrialsCommonBL{r,c,i}]),numTimePoints,numFreqPoints,numElecs),3)));
                dTFAcrossTrialsLogElectrodesCommonBL{r,c,i} = squeeze(nanmean(reshape(cell2mat([dTFAcrossTrialsCommonBL{r,c,i}]),numTimePoints,numFreqPoints,numElecs),3));
                dTFAcrossLogTrialsElectrodesCommonBL{r,c,i} = squeeze(nanmean(reshape(cell2mat([dTFAcrossLogTrialsCommonBL{r,c,i}]),numTimePoints,numFreqPoints,numElecs),3));
            end

            %____________________TF Measures_______________________
            fPos = params.freqVals>=params.fBandLow & params.freqVals<=params.fBandHigh; % frequency indices for the band of interest
            if isfield(params,'extraFreqBands')
                for l=1:length(params.extraFreqBands)
                    fPosExtra{l} = params.freqVals>=params.extraFreqBands{l}(1) & params.freqVals<=params.extraFreqBands{l}(2);
                end
            end
            
%             % Power in the selected frequency band
%             % Time course of the difference power in the selected band
%             dBandPowerAcrossTrials{r,c,i} = cellfun(@(x) sum(x(:,fPos),2),dTFAcrossTrials{r,c,i},'UniformOutput',0);
%             dBandPowerAcrossLogTrials{r,c,i} = cellfun(@(x) sum(x(:,fPos),2),dTFAcrossLogTrials{r,c,i},'UniformOutput',0);
%             
%             dBandPowerAcrossTrialsCommonBL{r,c,i} = cellfun(@(x) sum(x(:,fPos),2),dTFAcrossTrialsCommonBL{r,c,i},'UniformOutput',0);
%             dBandPowerAcrossLogTrialsCommonBL{r,c,i} = cellfun(@(x) sum(x(:,fPos),2),dTFAcrossLogTrialsCommonBL{r,c,i},'UniformOutput',0);
%             
%             % Averaged across electrodes
%             dBandPowerAcrossTrialsElectrodes{r,c,i} = sum(dTFAcrossTrialsElectrodes{r,c,i}(:,fPos),2);
%             dBandPowerAcrossTrialsLogElectrodes{r,c,i} = sum(dTFAcrossTrialsLogElectrodes{r,c,i}(:,fPos),2);
%             dBandPowerAcrossLogTrialsElectrodes{r,c,i} = sum(dTFAcrossLogTrialsElectrodes{r,c,i}(:,fPos),2);
%             
%             dBandPowerAcrossTrialsElectrodesCommonBL{r,c,i} = sum(dTFAcrossTrialsElectrodesCommonBL{r,c,i}(:,fPos),2);
%             dBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,i} = sum(dTFAcrossTrialsLogElectrodesCommonBL{r,c,i}(:,fPos),2);
%             dBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,i} = sum(dTFAcrossLogTrialsElectrodesCommonBL{r,c,i}(:,fPos),2);
            
            % Power in the selected frequency band
            % Time course of the difference power in the selected band
            dBandPowerAcrossTrials{r,c,i} = cellfun(@(x,y) conv2Log(sum(x(:,fPos),2)./sum(y(:,fPos),2)),TFAcrossTrials{r,c,i},baselineTFAcrossTrials,'UniformOutput',0);
            dBandPowerAcrossLogTrials{r,c,i} = cellfun(@(x,y) sum(x(:,fPos),2) - sum(y(:,fPos),2),TFAcrossLogTrials{r,c,i},baselineTFAcrossLogTrials,'UniformOutput',0);
            
            dBandPowerAcrossTrialsCommonBL{r,c,i} = cellfun(@(x,y) conv2Log(sum(x(:,fPos),2)./sum(y(:,fPos),2)),TFAcrossTrials{r,c,i},commonBLAcrossTrials','UniformOutput',0);
            dBandPowerAcrossLogTrialsCommonBL{r,c,i} = cellfun(@(x,y) sum(x(:,fPos),2) - sum(y(:,fPos),2),TFAcrossLogTrials{r,c,i},commonBLAcrossLogTrials','UniformOutput',0);
            
            % Averaged across electrodes
            dBandPowerAcrossTrialsElectrodes{r,c,i} = conv2Log(nanmean(10.^cell2mat(dBandPowerAcrossTrials{r,c,i}),2));
            dBandPowerAcrossTrialsLogElectrodes{r,c,i} = nanmean(cell2mat(dBandPowerAcrossTrials{r,c,i}),2);
            dBandPowerAcrossLogTrialsElectrodes{r,c,i} = nanmean(cell2mat(dBandPowerAcrossLogTrials{r,c,i}),2);
            
            dBandPowerAcrossTrialsElectrodesCommonBL{r,c,i} = conv2Log(nanmean(10.^cell2mat(dBandPowerAcrossTrialsCommonBL{r,c,i}),2));
            dBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,i} = nanmean(cell2mat(dBandPowerAcrossTrialsCommonBL{r,c,i}),2);
            dBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,i} = nanmean(cell2mat(dBandPowerAcrossLogTrialsCommonBL{r,c,i}),2);
            
            % Peak frequency in the selected frequency band
            bandFreqVals = params.freqVals(fPos);
            
            clear peakPower peakPowerLog peakPowerCommonBL peakPowerLogCommonBL
            peakPower = cellfun(@(x) max(x(:,fPos),[],2), dTFAcrossTrials{r,c,i}, 'uniformoutput',0);
            peakPowerLog = cellfun(@(x) max(x(:,fPos),[],2), dTFAcrossLogTrials{r,c,i}, 'uniformoutput',0);
            peakPowerCommonBL = cellfun(@(x) max(x(:,fPos),[],2), dTFAcrossTrialsCommonBL{r,c,i}, 'uniformoutput',0);
            peakPowerLogCommonBL = cellfun(@(x) max(x(:,fPos),[],2), dTFAcrossLogTrialsCommonBL{r,c,i}, 'uniformoutput',0);
            for tx=1:numTimePoints
                peakFreqdBandPowerAcrossTrials{r,c,i}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPos)==y(tx)), dTFAcrossTrials{r,c,i}, peakPower, 'uniformoutput',0));
                peakFreqdBandPowerAcrossLogTrials{r,c,i}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPos)==y(tx)), dTFAcrossLogTrials{r,c,i}, peakPowerLog, 'uniformoutput',0));
                
                peakFreqdBandPowerAcrossTrialsCommonBL{r,c,i}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPos)==y(tx)), dTFAcrossTrialsCommonBL{r,c,i}, peakPowerCommonBL, 'uniformoutput',0));
                peakFreqdBandPowerAcrossLogTrialsCommonBL{r,c,i}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPos)==y(tx)), dTFAcrossLogTrialsCommonBL{r,c,i}, peakPowerLogCommonBL, 'uniformoutput',0));
            end
            
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common')
                clear weights
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:))); % number of repeats per electrode for a session
                peakFreqdBandPowerAcrossTrialsElectrodes{r,c,i} = peakFreqdBandPowerAcrossTrials{r,c,i}*weights./sum(weights);
                peakFreqdBandPowerAcrossLogTrialsElectrodes{r,c,i} = peakFreqdBandPowerAcrossLogTrials{r,c,i}*weights./sum(weights);

                peakFreqdBandPowerAcrossTrialsElectrodesCommonBL{r,c,i} = peakFreqdBandPowerAcrossTrialsCommonBL{r,c,i}*weights./sum(weights);
                peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,i} = peakFreqdBandPowerAcrossLogTrialsCommonBL{r,c,i}*weights./sum(weights);
            else
                peakFreqdBandPowerAcrossTrialsElectrodes{r,c,i} = mean(peakFreqdBandPowerAcrossTrials{r,c,i},2);
                peakFreqdBandPowerAcrossLogTrialsElectrodes{r,c,i} = mean(peakFreqdBandPowerAcrossLogTrials{r,c,i},2);

                peakFreqdBandPowerAcrossTrialsElectrodesCommonBL{r,c,i} = mean(peakFreqdBandPowerAcrossTrialsCommonBL{r,c,i},2);
                peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,i} = mean(peakFreqdBandPowerAcrossLogTrialsCommonBL{r,c,i},2);
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

    dBandPowerAcrossTrialsEB = cell(numRows,numCols,numProtocls,numBands);
    dBandPowerAcrossLogTrialsEB = cell(numRows,numCols,numProtocls,numBands);
    dBandPowerAcrossTrialsCommonBLEB = cell(numRows,numCols,numProtocls,numBands);
    dBandPowerAcrossLogTrialsCommonBLEB = cell(numRows,numCols,numProtocls,numBands);

    dBandPowerAcrossTrialsElectrodesEB = cell(numRows,numCols,numProtocls,numBands);
    dBandPowerAcrossTrialsLogElectrodesEB = cell(numRows,numCols,numProtocls,numBands);
    dBandPowerAcrossLogTrialsElectrodesEB = cell(numRows,numCols,numProtocls,numBands);

    dBandPowerAcrossTrialsElectrodesCommonBLEB = cell(numRows,numCols,numProtocls,numBands);
    dBandPowerAcrossTrialsLogElectrodesCommonBLEB = cell(numRows,numCols,numProtocls,numBands);
    dBandPowerAcrossLogTrialsElectrodesCommonBLEB = cell(numRows,numCols,numProtocls,numBands);

    peakFreqdBandPowerAcrossTrialsEB = cell(numRows,numCols,numProtocls,numBands);
    peakFreqdBandPowerAcrossLogTrialsEB = cell(numRows,numCols,numProtocls,numBands);

    peakFreqdBandPowerAcrossTrialsCommonBLEB = cell(numRows,numCols,numProtocls,numBands);
    peakFreqdBandPowerAcrossLogTrialsCommonBLEB = cell(numRows,numCols,numProtocls,numBands);

    peakFreqdBandPowerAcrossTrialsElectrodesEB = cell(numRows,numCols,numProtocls,numBands);
    peakFreqdBandPowerAcrossLogTrialsElectrodesEB = cell(numRows,numCols,numProtocls,numBands);

    peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB = cell(numRows,numCols,numProtocls,numBands);
    peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBLEB = cell(numRows,numCols,numProtocls,numBands);

    for r = 1:numRows
        for c = 1:numCols
            for i = 1:numProtocols
                for l=1:length(fPosExtra)
%                     % Power in the selected frequency band
%                     % Time course of the difference power in the selected band
%                     dBandPowerAcrossTrialsEB{r,c,i,l} = cellfun(@(x) sum(x(:,fPosExtra{l}),2),dTFAcrossTrials{r,c,i},'UniformOutput',0);
%                     dBandPowerAcrossLogTrialsEB{r,c,i,l} = cellfun(@(x) sum(x(:,fPosExtra{l}),2),dTFAcrossLogTrials{r,c,i},'UniformOutput',0);
% 
%                     dBandPowerAcrossTrialsCommonBLEB{r,c,i,l} = cellfun(@(x) sum(x(:,fPosExtra{l}),2),dTFAcrossTrialsCommonBL{r,c,i},'UniformOutput',0);
%                     dBandPowerAcrossLogTrialsCommonBLEB{r,c,i,l} = cellfun(@(x) sum(x(:,fPosExtra{l}),2),dTFAcrossLogTrialsCommonBL{r,c,i},'UniformOutput',0);
% 
%                     % Averaged across electrodes
%                     dBandPowerAcrossTrialsElectrodesEB{r,c,i,l} = sum(dTFAcrossTrialsElectrodes{r,c,i}(:,fPosExtra{l}),2);
%                     dBandPowerAcrossTrialsLogElectrodesEB{r,c,i,l} = sum(dTFAcrossTrialsLogElectrodes{r,c,i}(:,fPosExtra{l}),2);
%                     dBandPowerAcrossLogTrialsElectrodesEB{r,c,i,l} = sum(dTFAcrossLogTrialsElectrodes{r,c,i}(:,fPosExtra{l}),2);
% 
%                     dBandPowerAcrossTrialsElectrodesCommonBLEB{r,c,i,l} = sum(dTFAcrossTrialsElectrodesCommonBL{r,c,i}(:,fPosExtra{l}),2);
%                     dBandPowerAcrossTrialsLogElectrodesCommonBLEB{r,c,i,l} = sum(dTFAcrossTrialsLogElectrodesCommonBL{r,c,i}(:,fPosExtra{l}),2);
%                     dBandPowerAcrossLogTrialsElectrodesCommonBLEB{r,c,i,l} = sum(dTFAcrossLogTrialsElectrodesCommonBL{r,c,i}(:,fPosExtra{l}),2);

                    
                    % Power in the selected frequency band
                    % Time course of the difference power in the selected band
                    dBandPowerAcrossTrialsEB{r,c,i,l} = cellfun(@(x,y) conv2Log(sum(x(:,fPosExtra{l}),2)./sum(y(:,fPosExtra{l}),2)),TFAcrossTrials{r,c,i},baselineTFAcrossTrials,'UniformOutput',0);
                    dBandPowerAcrossLogTrialsEB{r,c,i,l} = cellfun(@(x,y) sum(x(:,fPosExtra{l}),2) - sum(y(:,fPosExtra{l}),2),TFAcrossLogTrials{r,c,i},baselineTFAcrossLogTrials,'UniformOutput',0);

                    dBandPowerAcrossTrialsCommonBLEB{r,c,i,l} = cellfun(@(x,y) conv2Log(sum(x(:,fPosExtra{l}),2)./sum(y(:,fPosExtra{l}),2)),TFAcrossTrials{r,c,i},commonBLAcrossTrials','UniformOutput',0);
                    dBandPowerAcrossLogTrialsCommonBLEB{r,c,i,l} = cellfun(@(x,y) sum(x(:,fPosExtra{l}),2) - sum(y(:,fPosExtra{l}),2),TFAcrossLogTrials{r,c,i},commonBLAcrossLogTrials','UniformOutput',0);
                    
                    % Averaged across electrodes
                    dBandPowerAcrossTrialsElectrodesEB{r,c,i,l} = conv2Log(nanmean(10.^cell2mat(dBandPowerAcrossTrialsEB{r,c,i,l}),2));
                    dBandPowerAcrossTrialsLogElectrodesEB{r,c,i,l} = nanmean(cell2mat(dBandPowerAcrossTrialsEB{r,c,i,l}),2);
                    dBandPowerAcrossLogTrialsElectrodesEB{r,c,i,l} = nanmean(cell2mat(dBandPowerAcrossLogTrialsEB{r,c,i,l}),2);

                    dBandPowerAcrossTrialsElectrodesCommonBLEB{r,c,i,l} = conv2Log(nanmean(10.^cell2mat(dBandPowerAcrossTrialsCommonBLEB{r,c,i,l}),2));
                    dBandPowerAcrossTrialsLogElectrodesCommonBLEB{r,c,i,l} = nanmean(cell2mat(dBandPowerAcrossTrialsCommonBLEB{r,c,i,l}),2);
                    dBandPowerAcrossLogTrialsElectrodesCommonBLEB{r,c,i,l} = nanmean(cell2mat(dBandPowerAcrossLogTrialsCommonBLEB{r,c,i,l}),2);
                    
                    % Peak frequency in the selected frequency band
                    bandFreqVals = params.freqVals(fPosExtra{l});

                    clear peakPower peakPowerLog peakPowerCommonBL peakPowerLogCommonBL
                    peakPower = cellfun(@(x) max(x(:,fPosExtra{l}),[],2), dTFAcrossTrials{r,c,i}, 'uniformoutput',0);
                    peakPowerLog = cellfun(@(x) max(x(:,fPosExtra{l}),[],2), dTFAcrossLogTrials{r,c,i}, 'uniformoutput',0);
                    peakPowerCommonBL = cellfun(@(x) max(x(:,fPosExtra{l}),[],2), dTFAcrossTrialsCommonBL{r,c,i}, 'uniformoutput',0);
                    peakPowerLogCommonBL = cellfun(@(x) max(x(:,fPosExtra{l}),[],2), dTFAcrossLogTrialsCommonBL{r,c,i}, 'uniformoutput',0);
                    for tx=1:numTimePoints
                        peakFreqdBandPowerAcrossTrialsEB{r,c,i,l}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPosExtra{l})==y(tx)), dTFAcrossTrials{r,c,i}, peakPower, 'uniformoutput',0));
                        peakFreqdBandPowerAcrossLogTrialsEB{r,c,i,l}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPosExtra{l})==y(tx)), dTFAcrossLogTrials{r,c,i}, peakPowerLog, 'uniformoutput',0));

                        peakFreqdBandPowerAcrossTrialsCommonBLEB{r,c,i,l}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPosExtra{l})==y(tx)), dTFAcrossTrialsCommonBL{r,c,i}, peakPowerCommonBL, 'uniformoutput',0));
                        peakFreqdBandPowerAcrossLogTrialsCommonBLEB{r,c,i,l}(tx,:) = cell2mat(cellfun(@(x,y) bandFreqVals(x(tx,fPosExtra{l})==y(tx)), dTFAcrossLogTrialsCommonBL{r,c,i}, peakPowerLogCommonBL, 'uniformoutput',0));
                    end

                    if params.weightedmean && ~strcmpi(params.badtrialsoption,'common')
                        clear weights
                        weights = cell2mat(squeeze(params.numRepeats(r,c,i,:))); % number of repeats per electrode for a session
                        peakFreqdBandPowerAcrossTrialsElectrodesEB{r,c,i,l} = peakFreqdBandPowerAcrossTrials{r,c,i}*weights./sum(weights);
                        peakFreqdBandPowerAcrossLogTrialsElectrodesEB{r,c,i,l} = peakFreqdBandPowerAcrossLogTrials{r,c,i}*weights./sum(weights);

                        peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB{r,c,i,l} = peakFreqdBandPowerAcrossTrialsCommonBL{r,c,i}*weights./sum(weights);
                        peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBLEB{r,c,i,l} = peakFreqdBandPowerAcrossLogTrialsCommonBL{r,c,i}*weights./sum(weights);
                    else
                        peakFreqdBandPowerAcrossTrialsElectrodesEB{r,c,i,l} = mean(peakFreqdBandPowerAcrossTrials{r,c,i},2);
                        peakFreqdBandPowerAcrossLogTrialsElectrodesEB{r,c,i,l} = mean(peakFreqdBandPowerAcrossLogTrials{r,c,i},2);

                        peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB{r,c,i,l} = mean(peakFreqdBandPowerAcrossTrials{r,c,i},2);
                        peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBLEB{r,c,i,l} = mean(peakFreqdBandPowerAcrossLogTrials{r,c,i},2);
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
% Check if plot handles have been passed, else create a new figure and plot
% handles if newfiguresflag is set
if ~isfield(params,'hTFFig') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hTFFig = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hTFFig') 
    % Raw TFs
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    meanTFAcrossTrialsElectrodesSessions(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(meanTFAcrossTrialsElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    meanTFAcrossTrialsLogElectrodesSessions(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(meanTFAcrossTrialsLogElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    meanTFAcrossLogTrialsElectrodesSessions(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(meanTFAcrossLogTrialsElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the total mean TF 
                    % spectrum
                else
                    meanTFAcrossTrialsElectrodesSessions(r,c,:,:) = squeeze(nanmean(cell2mat(meanTFAcrossTrialsElectrodes(r,c,:)),3));
                    meanTFAcrossTrialsLogElectrodesSessions(r,c,:,:) = squeeze(nanmean(cell2mat(meanTFAcrossTrialsLogElectrodes(r,c,:)),3));
                    meanTFAcrossLogTrialsElectrodesSessions(r,c,:,:) = squeeze(nanmean(cell2mat(meanTFAcrossLogTrialsElectrodes(r,c,:)),3));
                end
            else
                meanTFAcrossTrialsElectrodesSessions(r,c,:,:) = meanTFAcrossTrialsElectrodes{r,c,:};
                meanTFAcrossTrialsLogElectrodesSessions(r,c,:,:) = meanTFAcrossTrialsLogElectrodes{r,c,:};
                meanTFAcrossLogTrialsElectrodesSessions(r,c,:,:) = meanTFAcrossLogTrialsElectrodes{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hTFFig(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hTFFig(x,y));
            end
            pcolor(params.TFtimeVals,params.freqVals,squeeze(meanTFAcrossTrialsLogElectrodesSessions(r,c,:,:))');
            setAxesProperties(params.hTFFig(x,y),params,x,y,rows);
            showTitleAndN(params.hTFFig(x,y),params,r,c);
            shading interp; colormap(jet);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end
    set(params.hTFFig,'ylim',ylimits);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean TF across trials, electrodes and sessions','FitBoxToText','on');
    end
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hTFFig(r,c));
                delete(axisHandle);
            end
        end
    end 
end

%==========================================================================
% Difference TFs: change of power from baseline to stimulus period
if ~isfield(params,'hdiffTFFig') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hdiffTFFig = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hdiffTFFig')
    dTFAcrossTrialsElectrodesSessions = zeros(numRows,numCols,numTimePoints,numFreqPoints);
    dTFAcrossTrialsLogElectrodesSessions = zeros(numRows,numCols,numTimePoints,numFreqPoints);
    dTFAcrossLogTrialsElectrodesSessions = zeros(numRows,numCols,numTimePoints,numFreqPoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    dTFAcrossTrialsElectrodesSessions(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dTFAcrossTrialsElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    dTFAcrossTrialsLogElectrodesSessions(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dTFAcrossTrialsLogElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    dTFAcrossLogTrialsElectrodesSessions(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dTFAcrossLogTrialsElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the total dTF spectrum
                else
                    dTFAcrossTrialsElectrodesSessions(r,c,:,:) = squeeze(nanmean(cell2mat(dTFAcrossTrialsElectrodes(r,c,:)),3));
                    dTFAcrossTrialsLogElectrodesSessions(r,c,:,:) = squeeze(nanmean(cell2mat(dTFAcrossTrialsLogElectrodes(r,c,:)),3));
                    dTFAcrossLogTrialsElectrodesSessions(r,c,:,:) = squeeze(nanmean(cell2mat(dTFAcrossLogTrialsElectrodes(r,c,:)),3));
                end
            else
                dTFAcrossTrialsElectrodesSessions(r,c,:,:) = dTFAcrossTrialsElectrodes{r,c,:};
                dTFAcrossTrialsLogElectrodesSessions(r,c,:,:) = dTFAcrossTrialsLogElectrodes{r,c,:};
                dTFAcrossLogTrialsElectrodesSessions(r,c,:,:) = dTFAcrossLogTrialsElectrodes{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hdiffTFFig(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hdiffTFFig(x,y));
            end
            pcolor(params.TFtimeVals,params.freqVals,10*squeeze(dTFAcrossTrialsLogElectrodesSessions(r,c,:,:))');
            setAxesProperties(params.hdiffTFFig(x,y),params,x,y,rows);
            showTitleAndN(params.hdiffTFFig(x,y),params,r,c);
            shading interp; colormap(jet);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end
    set(params.hdiffTFFig,'clim',[params.cmin params.cmax]);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean TF across trials, electrodes and sessions','FitBoxToText','on');
    end
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hdiffTFFig(r,c));
                delete(axisHandle);
            end
        end
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Difference TFs: change of power from baseline to stimulus period
% With common baseline across conditions per electrode per session
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
    dTFAcrossTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints,numFreqPoints);
    dTFAcrossTrialsLogElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints,numFreqPoints);
    dTFAcrossLogTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints,numFreqPoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    dTFAcrossTrialsElectrodesSessionsCommonBL(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dTFAcrossTrialsElectrodesCommonBL(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dTFAcrossTrialsLogElectrodesCommonBL(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    dTFAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:,:) = squeeze(sum(reshape(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dTFAcrossLogTrialsElectrodesCommonBL(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),numTimePoints,numFreqPoints,numProtocols),3));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the total dTF spectrum
                else
                    dTFAcrossTrialsElectrodesSessionsCommonBL(r,c,:,:) = squeeze(nanmean(cell2mat(dTFAcrossTrialsElectrodesCommonBL(r,c,:)),3));
                    dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:) = squeeze(nanmean(cell2mat(dTFAcrossTrialsLogElectrodesCommonBL(r,c,:)),3));
                    dTFAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:,:) = squeeze(nanmean(cell2mat(dTFAcrossLogTrialsElectrodesCommonBL(r,c,:)),3));
                end
            else
                dTFAcrossTrialsElectrodesSessionsCommonBL(r,c,:,:) = dTFAcrossTrialsElectrodesCommonBL{r,c,:};
                dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:) = dTFAcrossTrialsLogElectrodesCommonBL{r,c,:};
                dTFAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:,:) = dTFAcrossLogTrialsElectrodesCommonBL{r,c,:};
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
            cLims = [params.cmin params.cmax];
            CData = 10*squeeze(dTFAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:,:))';
            spHandle = subplot(params.hdiffTFCommonBLFig(x,y)); 
            [imData] = plotPColorAsImage(spHandle,params.TFtimeVals,params.freqVals,CData,tfXLims,[],cLims);

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
if ~isfield(params,'hBandPowerWithTime') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hBandPowerWithTime = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hBandPowerWithTime')
    dBandPowerAcrossTrialsElectrodesSessions = zeros(numRows,numCols,numTimePoints);
    dBandPowerAcrossTrialsLogElectrodesSessions = zeros(numRows,numCols,numTimePoints);
    dBandPowerAcrossLogTrialsElectrodesSessions = zeros(numRows,numCols,numTimePoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    dBandPowerAcrossTrialsElectrodesSessions(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dBandPowerAcrossTrialsElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    dBandPowerAcrossTrialsLogElectrodesSessions(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dBandPowerAcrossTrialsLogElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    dBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dBandPowerAcrossLogTrialsElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the total dBandPower
                else
                    dBandPowerAcrossTrialsElectrodesSessions(r,c,:) = squeeze(nanmean(cell2mat(dBandPowerAcrossTrialsElectrodes(r,c,:)),3));
                    dBandPowerAcrossTrialsLogElectrodesSessions(r,c,:) = squeeze(nanmean(cell2mat(dBandPowerAcrossTrialsLogElectrodes(r,c,:)),3));
                    dBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = squeeze(nanmean(cell2mat(dBandPowerAcrossLogTrialsElectrodes(r,c,:)),3));
                end
            else
                dBandPowerAcrossTrialsElectrodesSessions(r,c,:) = dBandPowerAcrossTrialsElectrodes{r,c,:};
                dBandPowerAcrossTrialsLogElectrodesSessions(r,c,:) = dBandPowerAcrossTrialsLogElectrodes{r,c,:};
                dBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = dBandPowerAcrossLogTrialsElectrodes{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hBandPowerWithTime(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hBandPowerWithTime(x,y));
            end
            plot(params.TFtimeVals,10*squeeze(dBandPowerAcrossTrialsLogElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
            hold on;
            line([params.timeVals(1) params.timeVals(end)],[0 0],'color',[0.5 0.4 0.7]);
            setAxesProperties(params.hBandPowerWithTime(x,y),params,x,y,rows);
            showTitleAndN(params.hBandPowerWithTime(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end

    set(params.hBandPowerWithTime,'ylim',ylimits);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean TF across trials, electrodes and sessions','FitBoxToText','on');
    end
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hBandPowerWithTime(r,c));
                delete(axisHandle);
            end
        end
    end
end

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
    dBandPowerAcrossTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints);
    dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints);
    dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    dBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dBandPowerAcrossTrialsElectrodesCommonBL(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dBandPowerAcrossTrialsLogElectrodesCommonBL(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(dBandPowerAcrossLogTrialsElectrodesCommonBL(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the total dBandPower
                else
                    dBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(cell2mat(dBandPowerAcrossTrialsElectrodesCommonBL(r,c,:)),3));
                    dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(cell2mat(dBandPowerAcrossTrialsLogElectrodesCommonBL(r,c,:)),3));
                    dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(cell2mat(dBandPowerAcrossLogTrialsElectrodesCommonBL(r,c,:)),3));
                end
            else
                dBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = dBandPowerAcrossTrialsElectrodesCommonBL{r,c,:};
                dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = dBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,:};
                dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = dBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,:};
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
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
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
if ~isfield(params,'hPeakFreqWithTime') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hPeakFreqWithTime = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hPeakFreqWithTime')
    peakFreqdBandPowerAcrossTrialsElectrodesSessions = zeros(numRows,numCols,numTimePoints);
    peakFreqdBandPowerAcrossLogTrialsElectrodesSessions = zeros(numRows,numCols,numTimePoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    peakFreqdBandPowerAcrossTrialsElectrodesSessions(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(peakFreqdBandPowerAcrossTrialsElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    peakFreqdBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(peakFreqdBandPowerAcrossLogTrialsElectrodes(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the overall peak freq
                else
                    peakFreqdBandPowerAcrossTrialsElectrodesSessions(r,c,:) = squeeze(nanmean(cell2mat(peakFreqdBandPowerAcrossTrialsElectrodes(r,c,:)),3));
                    peakFreqdBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = squeeze(nanmean(cell2mat(peakFreqdBandPowerAcrossLogTrialsElectrodes(r,c,:)),3));
                end
            else
                peakFreqdBandPowerAcrossTrialsElectrodesSessions(r,c,:) = peakFreqdBandPowerAcrossTrialsElectrodes{r,c,:};
                peakFreqdBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = peakFreqdBandPowerAcrossLogTrialsElectrodes{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hPeakFreqWithTime(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hPeakFreqWithTime(x,y));
            end
            plot(params.TFtimeVals,squeeze(peakFreqdBandPowerAcrossTrialsElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
            hold on;
            line([params.timeVals(1) params.timeVals(end)],[0 0],'color',[0.5 0.4 0.7]);
            setAxesProperties(params.hPeakFreqWithTime(x,y),params,x,y,rows);
            showTitleAndN(params.hPeakFreqWithTime(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end

    set(params.hPeakFreqWithTime,'ylim',[params.fBandLow params.fBandHigh]);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean TF across trials, electrodes and sessions','FitBoxToText','on');
    end
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hPeakFreqWithTime(r,c));
                delete(axisHandle);
            end
        end
    end
end

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
    peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols,numTimePoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(peakFreqdBandPowerAcrossTrialsElectrodesCommonBL(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(sum(cell2mat(cellfun(@(x,y) x*y/sum(weights), reshape(peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL(r,c,:),1,[]), num2cell(weights'),'uniformoutput',0)),2));
                    % x*y/sum(weights) is the contribution from each session,
                    % and then we sum these to get the overall peak freq
                else
                    peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(cell2mat(peakFreqdBandPowerAcrossTrialsElectrodesCommonBL(r,c,:)),3));
                    peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(cell2mat(peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL(r,c,:)),3));
                end
            else
                peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = peakFreqdBandPowerAcrossTrialsElectrodesCommonBL{r,c,:};
                peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,:};
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
    if strcmpi(params.baselinestrategy,'conditionwise')
        switch params.tfstrategy
            case 'raw' % no logarithmic operation
                tfdata.TF = TFAcrossTrials;
                tfdata.TFMean = meanTFAcrossTrialsElectrodes;
                tfdata.diffTF = dTFAcrossTrials;
                tfdata.diffTFMean = dTFAcrossTrialsElectrodes;
                tfdata.dBandPower = dBandPowerAcrossTrials;
                tfdata.dBandPowerMean = dBandPowerAcrossTrialsElectrodes;
                tfdata.peakFreq = peakFreqdBandPowerAcrossTrials;
                tfdata.peakFreqMean = peakFreqdBandPowerAcrossTrialsElectrodes;
                if drawplots && newfiguresflag
                    tfdata.TFMeanSessions = meanTFAcrossTrialsElectrodesSessions;
                    tfdata.diffTFMeanSessions = dTFAcrossTrialsElectrodesSessions;
                    tfdata.dBandPowerMeanSessions = dBandPowerAcrossTrialsElectrodesSessions;
                    tfdata.peakFreqMeanSessions = peakFreqdBandPowerAcrossTrialsElectrodesSessions;
                end
                if isfield(params,'extraFreqBands')
                    tfdata.dBandPowerEB = dBandPowerAcrossTrialsEB;
                    tfdata.dBandPowerMeanEB = dBandPowerAcrossTrialsElectrodesEB;
                    tfdata.peakFreqEB = peakFreqdBandPowerAcrossTrialsEB;
                    tfdata.peakFreqMeanEB = peakFreqdBandPowerAcrossTrialsElectrodesEB;
                end

            case 'logTrial' % take log TF of every trial
                tfdata.TF = TFAcrossLogTrials;
                tfdata.TFMean = meanTFAcrossLogTrialsElectrodes;
                tfdata.diffTF = dTFAcrossLogTrials;
                tfdata.diffTFMean = dTFAcrossLogTrialsElectrodes;
                tfdata.dBandPower = dBandPowerAcrossLogTrials;
                tfdata.dBandPowerMean = dBandPowerAcrossLogTrialsElectrodes;
                tfdata.peakFreq = peakFreqdBandPowerAcrossLogTrials;
                tfdata.peakFreqMean = peakFreqdBandPowerAcrossLogTrialsElectrodes;
                if drawplots && newfiguresflag
                    tfdata.TFMeanSessions = meanTFAcrossLogTrialsElectrodesSessions;
                    tfdata.diffTFMeanSessions = dTFAcrossLogTrialsElectrodesSessions;
                    tfdata.dBandPowerMeanSessions = dBandPowerAcrossLogTrialsElectrodesSessions;
                    tfdata.peakFreqMeanSessions = peakFreqdBandPowerAcrossLogTrialsElectrodesSessions;
                end
                if isfield(params,'extraFreqBands')
                    tfdata.dBandPowerEB = dBandPowerAcrossLogTrialsEB;
                    tfdata.dBandPowerMeanEB = dBandPowerAcrossLogTrialsElectrodesEB;
                    tfdata.peakFreqEB = peakFreqdBandPowerAcrossLogTrialsEB;
                    tfdata.peakFreqMeanEB = peakFreqdBandPowerAcrossLogTrialsElectrodesEB;
                end
            otherwise % default - take log at the level of every electrode
                tfdata.TF = TFAcrossTrials;
                tfdata.TFMean = meanTFAcrossTrialsLogElectrodes;
                tfdata.diffTF = dTFAcrossTrials;
                tfdata.diffTFMean = dTFAcrossTrialsLogElectrodes;
                tfdata.dBandPower = dBandPowerAcrossTrials;
                tfdata.dBandPowerMean = dBandPowerAcrossTrialsLogElectrodes;
                tfdata.peakFreq = peakFreqdBandPowerAcrossTrials;
                tfdata.peakFreqMean = peakFreqdBandPowerAcrossTrialsElectrodes;
                if drawplots && newfiguresflag
                    tfdata.TFMeanSessions = meanTFAcrossTrialsLogElectrodesSessions;
                    tfdata.diffTFMeanSessions = dTFAcrossTrialsLogElectrodesSessions;
                    tfdata.dBandPowerMeanSessions = dBandPowerAcrossTrialsLogElectrodesSessions;
                    tfdata.peakFreqMeanSessions = peakFreqdBandPowerAcrossTrialsElectrodesSessions;
                end
                if isfield(params,'extraFreqBands')
                    tfdata.dBandPowerEB = dBandPowerAcrossTrialsEB;
                    tfdata.dBandPowerMeanEB = dBandPowerAcrossTrialsLogElectrodesEB;
                    tfdata.peakFreqEB = peakFreqdBandPowerAcrossTrialsEB;
                    tfdata.peakFreqMeanEB = peakFreqdBandPowerAcrossTrialsElectrodesEB;
                end
        end
    elseif strcmpi(params.baselinestrategy,'common')
        switch params.tfstrategy
            case 'raw' % no logarithmic operation
                tfdata.TF = TFAcrossTrials;
                tfdata.TFMean = meanTFAcrossTrialsElectrodes;
                tfdata.diffTF = dTFAcrossTrialsCommonBL;
                tfdata.diffTFMean = dTFAcrossTrialsElectrodesCommonBL;
                tfdata.dBandPower = dBandPowerAcrossTrialsCommonBL;
                tfdata.dBandPowerMean = dBandPowerAcrossTrialsElectrodesCommonBL;
                tfdata.peakFreq = peakFreqdBandPowerAcrossTrialsCommonBL;
                tfdata.peakFreqMean = peakFreqdBandPowerAcrossTrialsElectrodesCommonBL;
                if drawplots && newfiguresflag
                    tfdata.TFMeanSessions = meanTFAcrossTrialsElectrodesSessions;
                    tfdata.diffTFMeanSessions = dTFAcrossTrialsElectrodesSessionsCommonBL;
                    tfdata.dBandPowerMeanSessions = dBandPowerAcrossTrialsElectrodesSessionsCommonBL;
                    tfdata.peakFreqMeanSessions = peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL;
                end
                if isfield(params,'extraFreqBands')
                    tfdata.dBandPowerEB = dBandPowerAcrossTrialsCommonBLEB;
                    tfdata.dBandPowerMeanEB = dBandPowerAcrossTrialsElectrodesCommonBLEB;
                    tfdata.peakFreqEB = peakFreqdBandPowerAcrossTrialsCommonBLEB;
                    tfdata.peakFreqMeanEB = peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB;
                end

            case 'logTrial' % take log TF of every trial
                tfdata.TF = TFAcrossLogTrials;
                tfdata.TFMean = meanTFAcrossLogTrialsElectrodes;
                tfdata.diffTF = dTFAcrossLogTrialsCommonBL;
                tfdata.diffTFMean = dTFAcrossLogTrialsElectrodesCommonBL;
                tfdata.dBandPower = dBandPowerAcrossLogTrialsCommonBL;
                tfdata.dBandPowerMean = dBandPowerAcrossLogTrialsElectrodesCommonBL;
                tfdata.peakFreq = peakFreqdBandPowerAcrossLogTrialsCommonBL;
                tfdata.peakFreqMean = peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL;
                if drawplots && newfiguresflag
                    tfdata.TFMeanSessions = meanTFAcrossLogTrialsElectrodesSessions;
                    tfdata.diffTFMeanSessions = dTFAcrossLogTrialsElectrodesSessionsCommonBL;
                    tfdata.dBandPowerMeanSessions = dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL;
                    tfdata.peakFreqMeanSessions = peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL;
                end
                if isfield(params,'extraFreqBands')
                    tfdata.dBandPowerEB = dBandPowerAcrossLogTrialsCommonBLEB;
                    tfdata.dBandPowerMeanEB = dBandPowerAcrossLogTrialsElectrodesCommonBLEB;
                    tfdata.peakFreqEB = peakFreqdBandPowerAcrossLogTrialsCommonBLEB;
                    tfdata.peakFreqMeanEB = peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBLEB;
                end
            otherwise % default - take log at the level of every electrode
                tfdata.TF = TFAcrossTrials;
                tfdata.TFMean = meanTFAcrossTrialsLogElectrodes;
                tfdata.diffTF = dTFAcrossTrialsCommonBL;
                tfdata.diffTFMean = dTFAcrossTrialsLogElectrodesCommonBL;
                tfdata.dBandPower = dBandPowerAcrossTrialsCommonBL;
                tfdata.dBandPowerMean = dBandPowerAcrossTrialsLogElectrodesCommonBL;
                tfdata.peakFreq = peakFreqdBandPowerAcrossTrialsCommonBL;
                tfdata.peakFreqMean = peakFreqdBandPowerAcrossTrialsElectrodesCommonBL;
                if drawplots && newfiguresflag
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
        end
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