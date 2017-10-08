% plot mean Power Spectral Density (PSD) across electrodes, sessions, 
% conditions
% Vinay Shirhatti, 14 October 2016
% modified and developed from computeAndPlotMeanERPMeasures
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
%             the calculations and return the PSD data
% newfiguresflag: set this flag to 1 to plot everything. If the relevant 
%                 handles are passed then the plots are drawn in those
%                 subplots. Otherwise new figures are created and the plots
%                 are drawn in them.
% returnData : set this flag to 1 to return the relevant data, else psddata
%              is empty
%
% OUTPUTS
% params : modified params
% psddata : the calculated PSD data
% *************************************************************************

function [params,psddata] = computeAndPlotMeanPSDMeasures(data,params,drawplots,newfiguresflag,returnData)

if ~exist('drawplots','var'); drawplots=1; end
if ~exist('newfiguresflag','var'); newfiguresflag=0; end
if ~exist('returnData','var'); returnData=1; end

%--------------------------------------------------------------------------
% Read all sizes
numRows = size(data,1);
numCols = size(data,2);
numProtocols = size(data,3);
% numElecs = size(data{1},2);

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

%--------------------------------------------------------------------------
% Compute and plot PSDs
%--------------------------------------------------------------------------
% Initialize variables
PSDEveryTrialBL = cell(numRows,numCols,numProtocols);
PSDAcrossTrialsBL = cell(numRows,numCols,numProtocols);
serrPSDAcrossTrialsBL = cell(numRows,numCols,numProtocols);
meanPSDAcrossTrialsElectrodesBL = cell(numRows,numCols,numProtocols);

PSDEveryTrialST = cell(numRows,numCols,numProtocols);
PSDAcrossTrialsST = cell(numRows,numCols,numProtocols);
serrPSDAcrossTrialsST = cell(numRows,numCols,numProtocols);
meanPSDAcrossTrialsElectrodesST = cell(numRows,numCols,numProtocols);

PSDAcrossLogTrialsBL = cell(numRows,numCols,numProtocols);
PSDAcrossLogTrialsST = cell(numRows,numCols,numProtocols);

meanPSDAcrossTrialsLogElectrodesBL = cell(numRows,numCols,numProtocols);
meanPSDAcrossTrialsLogElectrodesST = cell(numRows,numCols,numProtocols);

meanPSDAcrossLogTrialsElectrodesBL = cell(numRows,numCols,numProtocols);
meanPSDAcrossLogTrialsElectrodesST = cell(numRows,numCols,numProtocols);

% get the baseline and stimulus indices
Fs = round(1/(params.timeVals(2)-params.timeVals(1)));
BLRange = uint16((params.blmax-params.blmin)*Fs);
blPos = find(params.timeVals>=params.blmin,1)+ (1:BLRange);
STRange = uint16((params.stmax-params.stmin)*Fs);
stPos = find(params.timeVals>=params.stmin,1)+ (1:STRange);

% Compute PSD and related measures
for i = 1:numProtocols
    for r = 1:numRows
        for c = 1:numCols
            disp([' PSD >> processing: protocol ' num2str(i) '/' num2str(numProtocols) ', row ' num2str(c) '/' num2str(numCols) ', col ' num2str(r) '/' num2str(numRows)]);
            
            clear dataBL dataST
            dataBL = cellfun(@(x) x(:,blPos),data{r,c,i},'UniformOutput',0); % baseline epoch
            dataST = cellfun(@(x) x(:,stPos),data{r,c,i},'UniformOutput',0); % stimulus epoch
            
            %[PSDEveryTrial{r,c,i},freqValsAll,serrPSDAcrossTrials{r,c,i}] = cellfun(@(y) mtspectrumc(y,params.mtmParams),cellfun(@(x) x',data{r,c,i},'UniformOutput',0),'UniformOutput',0);
            [PSDEveryTrialBL{r,c,i},~,serrPSDAcrossTrialsBL{r,c,i}] = cellfun(@(y) mtspectrumc(y,params.mtmParams),cellfun(@(x) x',dataBL,'UniformOutput',0),'UniformOutput',0);
            [PSDEveryTrialST{r,c,i},freqVals,serrPSDAcrossTrialsST{r,c,i}] = cellfun(@(y) mtspectrumc(y,params.mtmParams),cellfun(@(x) x',dataST,'UniformOutput',0),'UniformOutput',0);
            % cellfun(@(x) x',data{r,c,i},'UniformOutput',0) takes the
            % conditionwise signals for every electrode, every trial and
            % transposes them so that they are of the form: samples x
            % trials, as required by mtspectrumc. Thus the result of this
            % operation is a cell array (of size 1 x numElecs) consisting
            % of signals from all electrodes with each cell corresponding
            % to all trials from that electrode. mtspectrumc is then run on
            % each electrode's trials. Looks great! cellfun's lot of fun!
            % :D
            
            %PSDAcrossTrials{r,c,i} = cellfun(@(x) mean(x,2),PSDEveryTrial{r,c,i},'UniformOutput',0);
            PSDAcrossTrialsBL{r,c,i} = cellfun(@(x) mean(x,2),PSDEveryTrialBL{r,c,i},'UniformOutput',0);
            PSDAcrossTrialsST{r,c,i} = cellfun(@(x) mean(x,2),PSDEveryTrialST{r,c,i},'UniformOutput',0);
            % mean across trials for every electrode => mean PSD for each
            % electrode for this condition
            
            % Similarly mean across logarithm of every trial
            PSDAcrossLogTrialsBL{r,c,i} = cellfun(@(x) mean(conv2Log(x),2),PSDEveryTrialBL{r,c,i},'UniformOutput',0);
            PSDAcrossLogTrialsST{r,c,i} = cellfun(@(x) mean(conv2Log(x),2),PSDEveryTrialST{r,c,i},'UniformOutput',0);
            
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common')
                % weighted mean: the number of repeats per electrode may 
                % vary. In such a case one may take a weighted mean 
                % This approach is required only when you have unequal
                % repeats per electrode i.e. when the badTrials selection
                % strategy is not the 'common' badTrials strategy
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:)));
                meanPSDAcrossTrialsElectrodesBL{r,c,i} = cell2mat([PSDAcrossTrialsBL{r,c,i}])*weights/sum(weights); % this is a matrix multiplication
                meanPSDAcrossTrialsElectrodesST{r,c,i} = cell2mat([PSDAcrossTrialsST{r,c,i}])*weights/sum(weights);
                % weights dim: numElecs x 1 => contains the number of 
                % repeats for each electrode for this condition and session
                % cell2mat([PSDAcrossTrialsBL{r,c,i}]) dim: numFreqPoints x numElecs =>
                % contains the mean across trials for each electrode.
                % Through this matrix operation we weight each PSD by the
                % number of repeats used to compute it and then take the
                % mean of all PSDs across electrodes. The resultant matrix
                % is of dimensions numFreqPoints x 1 i.e. the weighted
                % average of PSDs across electrodes
                meanPSDAcrossTrialsLogElectrodesBL{r,c,i} = conv2Log(cell2mat([PSDAcrossTrialsBL{r,c,i}]))*weights/sum(weights); % this is a matrix multiplication
                meanPSDAcrossTrialsLogElectrodesST{r,c,i} = conv2Log(cell2mat([PSDAcrossTrialsST{r,c,i}]))*weights/sum(weights);
                
                meanPSDAcrossLogTrialsElectrodesBL{r,c,i} = cell2mat([PSDAcrossLogTrialsBL{r,c,i}])*weights/sum(weights); % this is a matrix multiplication
                meanPSDAcrossLogTrialsElectrodesST{r,c,i} = cell2mat([PSDAcrossLogTrialsST{r,c,i}])*weights/sum(weights);
            else
                %meanPSDAcrossTrialsElectrodes{r,c,i} = mean(cell2mat([PSDAcrossTrials{r,c,i}]),2);
                meanPSDAcrossTrialsElectrodesBL{r,c,i} = mean(cell2mat([PSDAcrossTrialsBL{r,c,i}]),2);
                meanPSDAcrossTrialsElectrodesST{r,c,i} = mean(cell2mat([PSDAcrossTrialsST{r,c,i}]),2);
                % mean across electrodes for each condition, each session

                % mean across logarithm of mean electrode PSDs
                meanPSDAcrossTrialsLogElectrodesBL{r,c,i} = mean(conv2Log(cell2mat([PSDAcrossTrialsBL{r,c,i}])),2);
                meanPSDAcrossTrialsLogElectrodesST{r,c,i} = mean(conv2Log(cell2mat([PSDAcrossTrialsST{r,c,i}])),2);

                % mean across mean electrode PSDs calculated from logarithm of every trial
                meanPSDAcrossLogTrialsElectrodesBL{r,c,i} = mean(cell2mat([PSDAcrossLogTrialsBL{r,c,i}]),2);
                meanPSDAcrossLogTrialsElectrodesST{r,c,i} = mean(cell2mat([PSDAcrossLogTrialsST{r,c,i}]),2);
            end
            
        end
    end
end

commonBaselinePSDAcrossTrials = cell(1,numProtocols);
commonBaselinePSDAcrossLogTrials = cell(1,numProtocols);
for i = 1:numProtocols % find the common baseline for each session separately
    % Compute common baseline across conditions per electrode
    numElecs = length(data{1,1,i});
    commonBaselinePSDAcrossTrials{i} = cell(1,numElecs);
    commonBaselinePSDAcrossLogTrials{i} = cell(1,numElecs);
    PSDAcrossTrialsBLTranspose = cellfun(@(x) x',PSDAcrossTrialsBL,'UniformOutput',0); % to align the electrodes along the rows
    PSDAcrossLogTrialsBLTranspose = cellfun(@(x) x',PSDAcrossLogTrialsBL,'UniformOutput',0); % to align the electrodes along the rows
    
    clear PSDAcrossTrialsBLAllConditions
    PSDAcrossTrialsBLAllConditions = [PSDAcrossTrialsBLTranspose{:,:,i}];
    PSDAcrossLogTrialsBLAllConditions = [PSDAcrossLogTrialsBLTranspose{:,:,i}];
    for en = 1:numElecs
        commonBaselinePSDAcrossTrials{i}{en} = mean(cell2mat(PSDAcrossTrialsBLAllConditions(en,:)),2);
        commonBaselinePSDAcrossLogTrials{i}{en} = mean(cell2mat(PSDAcrossLogTrialsBLAllConditions(en,:)),2);
    end
    % take all conditions together for a protocol and find the mean
    % baseline PSD across all conditions for each electrode separately.
    % This is done by summing the PSDs across conditions and dividing by
    % the total number of conditions for each electrode here
end

% Calculate PSD measures: difference PSD, band power
dPSDAcrossTrials = cell(numRows,numCols,numProtocols);
dPSDAcrossLogTrials = cell(numRows,numCols,numProtocols);

dPSDAcrossTrialsElectrodes = cell(numRows,numCols,numProtocols);
dPSDAcrossTrialsLogElectrodes = cell(numRows,numCols,numProtocols);
dPSDAcrossLogTrialsElectrodes = cell(numRows,numCols,numProtocols);

dPSDAcrossTrialsCommonBL = cell(numRows,numCols,numProtocols);
dPSDAcrossLogTrialsCommonBL = cell(numRows,numCols,numProtocols);

dPSDAcrossTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);
dPSDAcrossTrialsLogElectrodesCommonBL = cell(numRows,numCols,numProtocols);
dPSDAcrossLogTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);

dBandPowerAcrossTrials = cell(numRows,numCols,numProtocols);
dBandPowerAcrossLogTrials = cell(numRows,numCols,numProtocols);

dBandPowerAcrossTrialsElectrodes = cell(numRows,numCols,numProtocols);
dBandPowerAcrossTrialsLogElectrodes = cell(numRows,numCols,numProtocols);
dBandPowerAcrossLogTrialsElectrodes = cell(numRows,numCols,numProtocols);

dBandPowerAcrossTrialsCommonBL = cell(numRows,numCols,numProtocols);
dBandPowerAcrossLogTrialsCommonBL = cell(numRows,numCols,numProtocols);

dBandPowerAcrossTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);
dBandPowerAcrossTrialsLogElectrodesCommonBL = cell(numRows,numCols,numProtocols);
dBandPowerAcrossLogTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);

peakFreqdBandPowerAcrossTrials = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossLogTrials= cell(numRows,numCols,numProtocols);

peakFreqdBandPowerAcrossTrialsElectrodes = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossTrialsLogElectrodes = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossLogTrialsElectrodes = cell(numRows,numCols,numProtocols);

peakFreqdBandPowerAcrossTrialsCommonBL = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossLogTrialsCommonBL = cell(numRows,numCols,numProtocols);

peakFreqdBandPowerAcrossTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossTrialsLogElectrodesCommonBL = cell(numRows,numCols,numProtocols);
peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL = cell(numRows,numCols,numProtocols);

for r = 1:numRows
    for c = 1:numCols
        for i = 1:numProtocols
            %______________Difference PSDs (logarithmic)___________________
            % per electrode averaged across its trials
            dPSDAcrossTrials{r,c,i} = cellfun(@(x,y) conv2Log(x./y),PSDAcrossTrialsST{r,c,i},PSDAcrossTrialsBL{r,c,i},'UniformOutput',0);
            dPSDAcrossLogTrials{r,c,i} = cellfun(@(x,y) (x-y),PSDAcrossLogTrialsST{r,c,i},PSDAcrossLogTrialsBL{r,c,i},'UniformOutput',0);
            
            % calculate mean across electrodes; if weighted average is
            % opted for then the weighing is already considered while
            % calculating the meanPSDs => no further weighing required here
            dPSDAcrossTrialsElectrodes{r,c,i} = conv2Log(meanPSDAcrossTrialsElectrodesST{r,c,i}./meanPSDAcrossTrialsElectrodesBL{r,c,i});
            dPSDAcrossTrialsLogElectrodes{r,c,i} = meanPSDAcrossTrialsLogElectrodesST{r,c,i} - meanPSDAcrossTrialsLogElectrodesBL{r,c,i};
            dPSDAcrossLogTrialsElectrodes{r,c,i} = meanPSDAcrossLogTrialsElectrodesST{r,c,i} - meanPSDAcrossLogTrialsElectrodesBL{r,c,i};
            
            %__Difference PSDs(log) wrt common baseline across conditions__
            dPSDAcrossTrialsCommonBL{r,c,i} = cellfun(@(x,y) conv2Log(x./y), PSDAcrossTrialsST{r,c,i}',commonBaselinePSDAcrossTrials{i}','UniformOutput',0);
            dPSDAcrossLogTrialsCommonBL{r,c,i} = cellfun(@(x,y) x-y, PSDAcrossLogTrialsST{r,c,i}',commonBaselinePSDAcrossLogTrials{i}','UniformOutput',0);
            
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common') % averaging across electrodes with possibly unequal repeats
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:))); % number of repeats per electrode for a session
                dPSDAcrossTrialsElectrodesCommonBL{r,c,i} = conv2Log(10.^cell2mat(dPSDAcrossTrialsCommonBL{r,c,i}')*weights./sum(weights));
                dPSDAcrossTrialsLogElectrodesCommonBL{r,c,i} = cell2mat(dPSDAcrossTrialsCommonBL{r,c,i}')*weights./sum(weights);
                dPSDAcrossLogTrialsElectrodesCommonBL{r,c,i} = cell2mat(dPSDAcrossLogTrialsCommonBL{r,c,i}')*weights./sum(weights);
                % This is a matrix multiplication
                % cell2mat(dPSDAcrossTrialsCommonBL{r,c,i}'): numFreqVals x numElecs;
                % weights: numElecs x 1
            else
                dPSDAcrossTrialsElectrodesCommonBL{r,c,i} = conv2Log(mean(10.^cell2mat(dPSDAcrossTrialsCommonBL{r,c,i}'),2));
                dPSDAcrossTrialsLogElectrodesCommonBL{r,c,i} = mean(cell2mat(dPSDAcrossTrialsCommonBL{r,c,i}'),2);
                dPSDAcrossLogTrialsElectrodesCommonBL{r,c,i} = mean(cell2mat(dPSDAcrossLogTrialsCommonBL{r,c,i}'),2);
            end
            
            params.freqVals = freqVals{1}; % store the freqVals

            %____________________PSD Measures_______________________
            fPos = params.freqVals>=params.fBandLow & params.freqVals<=params.fBandHigh; % frequency indices for the band of interest
            
            % Power in the selected frequency band
%             dBandPowerAcrossTrials{r,c,i} = cellfun(@(x) sum(x(fPos)),dPSDAcrossTrials{r,c,i},'UniformOutput',0);
%             dBandPowerAcrossLogTrials{r,c,i} = cellfun(@(x) sum(x(fPos)),dPSDAcrossLogTrials{r,c,i},'UniformOutput',0);
%             
%             dBandPowerAcrossTrialsElectrodes{r,c,i} = sum(dPSDAcrossTrialsElectrodes{r,c,i}(fPos));
%             dBandPowerAcrossTrialsLogElectrodes{r,c,i} = sum(dPSDAcrossTrialsLogElectrodes{r,c,i}(fPos));
%             dBandPowerAcrossLogTrialsElectrodes{r,c,i} = sum(dPSDAcrossLogTrialsElectrodes{r,c,i}(fPos));
%             
%             dBandPowerAcrossTrialsCommonBL{r,c,i} = cellfun(@(x) sum(x(fPos)),dPSDAcrossTrialsCommonBL{r,c,i},'UniformOutput',0);
%             dBandPowerAcrossLogTrialsCommonBL{r,c,i} = cellfun(@(x) sum(x(fPos)),dPSDAcrossLogTrialsCommonBL{r,c,i},'UniformOutput',0);
%             
%             dBandPowerAcrossTrialsElectrodesCommonBL{r,c,i} = sum(dPSDAcrossTrialsElectrodesCommonBL{r,c,i}(fPos));
%             dBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,i} = sum(dPSDAcrossTrialsLogElectrodesCommonBL{r,c,i}(fPos));
%             dBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,i} = sum(dPSDAcrossLogTrialsElectrodesCommonBL{r,c,i}(fPos));
            
            dBandPowerAcrossTrials{r,c,i} = cellfun(@(x,y) conv2Log(sum(x(fPos))./sum(y(fPos))),PSDAcrossTrialsST{r,c,i},PSDAcrossTrialsBL{r,c,i},'UniformOutput',0);
            dBandPowerAcrossLogTrials{r,c,i} = cellfun(@(x,y) sum(x(fPos))-sum(y(fPos)),PSDAcrossLogTrialsST{r,c,i},PSDAcrossLogTrialsBL{r,c,i},'UniformOutput',0);
            
            dBandPowerAcrossTrialsElectrodes{r,c,i} = conv2Log(nanmean(10.^cell2mat(dBandPowerAcrossTrials{r,c,i})));
            dBandPowerAcrossTrialsLogElectrodes{r,c,i} = nanmean(cell2mat(dBandPowerAcrossTrials{r,c,i}));
            dBandPowerAcrossLogTrialsElectrodes{r,c,i} = nanmean(cell2mat(dBandPowerAcrossLogTrials{r,c,i}));
            
            dBandPowerAcrossTrialsCommonBL{r,c,i} = cellfun(@(x,y) conv2Log(sum(x(fPos))./sum(y(fPos))),PSDAcrossTrialsST{r,c,i}',commonBaselinePSDAcrossTrials{i}','UniformOutput',0);
            dBandPowerAcrossLogTrialsCommonBL{r,c,i} = cellfun(@(x,y) sum(x(fPos))-sum(y(fPos)),PSDAcrossLogTrialsST{r,c,i}',commonBaselinePSDAcrossLogTrials{i}','UniformOutput',0);
            
            dBandPowerAcrossTrialsElectrodesCommonBL{r,c,i} = conv2Log(nanmean(10.^cell2mat(dBandPowerAcrossTrialsCommonBL{r,c,i})));
            dBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,i} = nanmean(cell2mat(dBandPowerAcrossTrialsCommonBL{r,c,i}));
            dBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,i} = nanmean(cell2mat(dBandPowerAcrossLogTrialsCommonBL{r,c,i}));
            
            % Peak frequency in the selected frequency band
            bandFreqVals = params.freqVals(fPos);
            
            peakFreqdBandPowerAcrossTrials{r,c,i} = cellfun(@(x) bandFreqVals(x(fPos)==max(x(fPos))), dPSDAcrossTrials{r,c,i},'UniformOutput',0);
            peakFreqdBandPowerAcrossLogTrials{r,c,i} = cellfun(@(x) bandFreqVals(x(fPos)==max(x(fPos))), dPSDAcrossLogTrials{r,c,i},'UniformOutput',0);
            
            peakFreqdBandPowerAcrossTrialsElectrodes{r,c,i} = bandFreqVals(dPSDAcrossTrialsElectrodes{r,c,i}(fPos)==max(dPSDAcrossTrialsElectrodes{r,c,i}(fPos)));
            peakFreqdBandPowerAcrossTrialsLogElectrodes{r,c,i} = bandFreqVals(dPSDAcrossTrialsLogElectrodes{r,c,i}(fPos)==max(dPSDAcrossTrialsLogElectrodes{r,c,i}(fPos)));
            peakFreqdBandPowerAcrossLogTrialsElectrodes{r,c,i} = bandFreqVals(dPSDAcrossLogTrialsElectrodes{r,c,i}(fPos)==max(dPSDAcrossLogTrialsElectrodes{r,c,i}(fPos)));
            
            peakFreqdBandPowerAcrossTrialsCommonBL{r,c,i} = cellfun(@(x) bandFreqVals(x(fPos)==max(x(fPos))), dPSDAcrossTrialsCommonBL{r,c,i},'UniformOutput',0);
            peakFreqdBandPowerAcrossLogTrialsCommonBL{r,c,i} = cellfun(@(x) bandFreqVals(x(fPos)==max(x(fPos))), dPSDAcrossLogTrialsCommonBL{r,c,i},'UniformOutput',0);
            
            peakFreqdBandPowerAcrossTrialsElectrodesCommonBL{r,c,i} = bandFreqVals(dPSDAcrossTrialsElectrodesCommonBL{r,c,i}(fPos)==max(dPSDAcrossTrialsElectrodesCommonBL{r,c,i}(fPos)));
            peakFreqdBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,i} = bandFreqVals(dPSDAcrossTrialsLogElectrodesCommonBL{r,c,i}(fPos)==max(dPSDAcrossTrialsLogElectrodesCommonBL{r,c,i}(fPos)));
            peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,i} = bandFreqVals(dPSDAcrossLogTrialsElectrodesCommonBL{r,c,i}(fPos)==max(dPSDAcrossLogTrialsElectrodesCommonBL{r,c,i}(fPos)));
            
        end
    end
end
numFreqPoints = length(params.freqVals);


if isfield(params,'extraFreqBands')

    for l=1:length(params.extraFreqBands)
        fPosExtra{l} = params.freqVals>=params.extraFreqBands{l}(1) & params.freqVals<=params.extraFreqBands{l}(2);
    end
    
    numBands = length(params.extraFreqBands);

    dBandPowerAcrossTrialsEB = cell(numRows,numCols,numProtocols,numBands);
    dBandPowerAcrossLogTrialsEB = cell(numRows,numCols,numProtocols,numBands);

    dBandPowerAcrossTrialsElectrodesEB = cell(numRows,numCols,numProtocols,numBands);
    dBandPowerAcrossTrialsLogElectrodesEB = cell(numRows,numCols,numProtocols,numBands);
    dBandPowerAcrossLogTrialsElectrodesEB = cell(numRows,numCols,numProtocols,numBands);

    dBandPowerAcrossTrialsCommonBLEB = cell(numRows,numCols,numProtocols,numBands);
    dBandPowerAcrossLogTrialsCommonBLEB = cell(numRows,numCols,numProtocols,numBands);

    dBandPowerAcrossTrialsElectrodesCommonBLEB = cell(numRows,numCols,numProtocols,numBands);
    dBandPowerAcrossTrialsLogElectrodesCommonBLEB = cell(numRows,numCols,numProtocols,numBands);
    dBandPowerAcrossLogTrialsElectrodesCommonBLEB = cell(numRows,numCols,numProtocols,numBands);

    peakFreqdBandPowerAcrossTrialsEB = cell(numRows,numCols,numProtocols,numBands);
    peakFreqdBandPowerAcrossLogTrialsEB = cell(numRows,numCols,numProtocols,numBands);

    peakFreqdBandPowerAcrossTrialsElectrodesEB = cell(numRows,numCols,numProtocols,numBands);
    peakFreqdBandPowerAcrossTrialsLogElectrodesEB = cell(numRows,numCols,numProtocols,numBands);
    peakFreqdBandPowerAcrossLogTrialsElectrodesEB = cell(numRows,numCols,numProtocols,numBands);

    peakFreqdBandPowerAcrossTrialsCommonBLEB = cell(numRows,numCols,numProtocols,numBands);
    peakFreqdBandPowerAcrossLogTrialsCommonBLEB = cell(numRows,numCols,numProtocols,numBands);

    peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB = cell(numRows,numCols,numProtocols,numBands);
    peakFreqdBandPowerAcrossTrialsLogElectrodesCommonBLEB = cell(numRows,numCols,numProtocols,numBands);
    peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBLEB = cell(numRows,numCols,numProtocols,numBands);
    
    for r = 1:numRows
        for c = 1:numCols
            for i = 1:numProtocols
                for l=1:length(fPosExtra)
                    % Power in the selected frequency band
                    dBandPowerAcrossTrialsEB{r,c,i,l} = cellfun(@(x,y) conv2Log(sum(x(fPosExtra{l}))./sum(y(fPosExtra{l}))),PSDAcrossTrialsST{r,c,i},PSDAcrossTrialsBL{r,c,i},'UniformOutput',0);
                    dBandPowerAcrossLogTrialsEB{r,c,i,l} = cellfun(@(x,y) sum(x(fPosExtra{l}))-sum(y(fPosExtra{l})),PSDAcrossLogTrialsST{r,c,i},PSDAcrossLogTrialsBL{r,c,i},'UniformOutput',0);

                    dBandPowerAcrossTrialsElectrodesEB{r,c,i,l} = conv2Log(nanmean(10.^cell2mat(dBandPowerAcrossTrialsEB{r,c,i,l})));
                    dBandPowerAcrossTrialsLogElectrodesEB{r,c,i,l} = nanmean(cell2mat(dBandPowerAcrossTrialsEB{r,c,i,l}));
                    dBandPowerAcrossLogTrialsElectrodesEB{r,c,i,l} = nanmean(cell2mat(dBandPowerAcrossLogTrialsEB{r,c,i,l}));

                    dBandPowerAcrossTrialsCommonBLEB{r,c,i,l} = cellfun(@(x,y) conv2Log(sum(x(fPosExtra{l}))./sum(y(fPosExtra{l}))),PSDAcrossTrialsST{r,c,i}',commonBaselinePSDAcrossTrials{i}','UniformOutput',0);
                    dBandPowerAcrossLogTrialsCommonBLEB{r,c,i,l} = cellfun(@(x,y) sum(x(fPosExtra{l}))-sum(y(fPosExtra{l})),PSDAcrossLogTrialsST{r,c,i}',commonBaselinePSDAcrossLogTrials{i}','UniformOutput',0);

                    dBandPowerAcrossTrialsElectrodesCommonBLEB{r,c,i,l} = conv2Log(nanmean(10.^cell2mat(dBandPowerAcrossTrialsCommonBLEB{r,c,i,l})));
                    dBandPowerAcrossTrialsLogElectrodesCommonBLEB{r,c,i,l} = nanmean(cell2mat(dBandPowerAcrossTrialsCommonBLEB{r,c,i,l}));
                    dBandPowerAcrossLogTrialsElectrodesCommonBLEB{r,c,i,l} = nanmean(cell2mat(dBandPowerAcrossLogTrialsCommonBLEB{r,c,i,l}));

                    % Peak frequency in the selected frequency band
                    bandFreqVals = params.freqVals(fPosExtra{l});

                    peakFreqdBandPowerAcrossTrialsEB{r,c,i,l} = cellfun(@(x) bandFreqVals(x(fPosExtra{l})==max(x(fPosExtra{l}))), dPSDAcrossTrials{r,c,i},'UniformOutput',0);
                    peakFreqdBandPowerAcrossLogTrialsEB{r,c,i,l} = cellfun(@(x) bandFreqVals(x(fPosExtra{l})==max(x(fPosExtra{l}))), dPSDAcrossLogTrials{r,c,i},'UniformOutput',0);

                    peakFreqdBandPowerAcrossTrialsElectrodesEB{r,c,i,l} = bandFreqVals(dPSDAcrossTrialsElectrodes{r,c,i}(fPosExtra{l})==max(dPSDAcrossTrialsElectrodes{r,c,i}(fPosExtra{l})));
                    peakFreqdBandPowerAcrossTrialsLogElectrodesEB{r,c,i,l} = bandFreqVals(dPSDAcrossTrialsLogElectrodes{r,c,i}(fPosExtra{l})==max(dPSDAcrossTrialsLogElectrodes{r,c,i}(fPosExtra{l})));
                    peakFreqdBandPowerAcrossLogTrialsElectrodesEB{r,c,i,l} = bandFreqVals(dPSDAcrossLogTrialsElectrodes{r,c,i}(fPosExtra{l})==max(dPSDAcrossLogTrialsElectrodes{r,c,i}(fPosExtra{l})));

                    peakFreqdBandPowerAcrossTrialsCommonBLEB{r,c,i,l} = cellfun(@(x) bandFreqVals(x(fPosExtra{l})==max(x(fPosExtra{l}))), dPSDAcrossTrialsCommonBL{r,c,i},'UniformOutput',0);
                    peakFreqdBandPowerAcrossLogTrialsCommonBLEB{r,c,i,l} = cellfun(@(x) bandFreqVals(x(fPosExtra{l})==max(x(fPosExtra{l}))), dPSDAcrossLogTrialsCommonBL{r,c,i},'UniformOutput',0);

                    peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB{r,c,i,l} = bandFreqVals(dPSDAcrossTrialsElectrodesCommonBL{r,c,i}(fPosExtra{l})==max(dPSDAcrossTrialsElectrodesCommonBL{r,c,i}(fPosExtra{l})));
                    peakFreqdBandPowerAcrossTrialsLogElectrodesCommonBLEB{r,c,i,l} = bandFreqVals(dPSDAcrossTrialsLogElectrodesCommonBL{r,c,i}(fPosExtra{l})==max(dPSDAcrossTrialsLogElectrodesCommonBL{r,c,i}(fPosExtra{l})));
                    peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBLEB{r,c,i,l} = bandFreqVals(dPSDAcrossLogTrialsElectrodesCommonBL{r,c,i}(fPosExtra{l})==max(dPSDAcrossLogTrialsElectrodesCommonBL{r,c,i}(fPosExtra{l})));
                end
            end
        end
    end
end

%==========================================================================
if drawplots
%==========================================================================
%--------------------------------------------------------------------------
% Draw the plots
%--------------------------------------------------------------------------
% Plot PSDs
% Raw PSDs: plot baseline and stimulus period PSDs
% Check if plot handles have been passed, else create a new figure and plot
% handles
if ~isfield(params,'hPSDFig') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hPSDFig = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hPSDFig')
    meanPSDAcrossTrialsElectrodesSessionsBL = zeros(numRows,numCols,numFreqPoints);
    meanPSDAcrossTrialsElectrodesSessionsST = zeros(numRows,numCols,numFreqPoints);
    meanPSDAcrossTrialsLogElectrodesSessionsBL = zeros(numRows,numCols,numFreqPoints);
    meanPSDAcrossTrialsLogElectrodesSessionsST = zeros(numRows,numCols,numFreqPoints);
    meanPSDAcrossLogTrialsElectrodesSessionsBL = zeros(numRows,numCols,numFreqPoints);
    meanPSDAcrossLogTrialsElectrodesSessionsST = zeros(numRows,numCols,numFreqPoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';

                    meanPSDAcrossTrialsElectrodesSessionsBL(r,c,:) = squeeze((squeeze(cell2mat(meanPSDAcrossTrialsElectrodesBL(r,c,:)))*weights)/sum(weights));
                    meanPSDAcrossTrialsElectrodesSessionsST(r,c,:) = squeeze((squeeze(cell2mat(meanPSDAcrossTrialsElectrodesST(r,c,:)))*weights)/sum(weights));

                    meanPSDAcrossTrialsLogElectrodesSessionsBL(r,c,:) = squeeze((squeeze(cell2mat(meanPSDAcrossTrialsLogElectrodesBL(r,c,:)))*weights)/sum(weights));
                    meanPSDAcrossTrialsLogElectrodesSessionsST(r,c,:) = squeeze((squeeze(cell2mat(meanPSDAcrossTrialsLogElectrodesST(r,c,:)))*weights)/sum(weights));

                    meanPSDAcrossLogTrialsElectrodesSessionsBL(r,c,:) = squeeze((squeeze(cell2mat(meanPSDAcrossLogTrialsElectrodesBL(r,c,:)))*weights)/sum(weights));
                    meanPSDAcrossLogTrialsElectrodesSessionsST(r,c,:) = squeeze((squeeze(cell2mat(meanPSDAcrossLogTrialsElectrodesST(r,c,:)))*weights)/sum(weights));
                else
                    meanPSDAcrossTrialsElectrodesSessionsBL(r,c,:) = squeeze(nanmean(squeeze(cell2mat(meanPSDAcrossTrialsElectrodesBL(r,c,:))),2));
                    meanPSDAcrossTrialsElectrodesSessionsST(r,c,:) = squeeze(nanmean(squeeze(cell2mat(meanPSDAcrossTrialsElectrodesST(r,c,:))),2));

                    meanPSDAcrossTrialsLogElectrodesSessionsBL(r,c,:) = squeeze(nanmean(squeeze(cell2mat(meanPSDAcrossTrialsLogElectrodesBL(r,c,:))),2));
                    meanPSDAcrossTrialsLogElectrodesSessionsST(r,c,:) = squeeze(nanmean(squeeze(cell2mat(meanPSDAcrossTrialsLogElectrodesST(r,c,:))),2));

                    meanPSDAcrossLogTrialsElectrodesSessionsBL(r,c,:) = squeeze(nanmean(squeeze(cell2mat(meanPSDAcrossLogTrialsElectrodesBL(r,c,:))),2));
                    meanPSDAcrossLogTrialsElectrodesSessionsST(r,c,:) = squeeze(nanmean(squeeze(cell2mat(meanPSDAcrossLogTrialsElectrodesST(r,c,:))),2));
                end
            else
                meanPSDAcrossTrialsElectrodesSessionsBL(r,c,:) = meanPSDAcrossTrialsElectrodesBL{r,c,:};
                meanPSDAcrossTrialsElectrodesSessionsST(r,c,:) = meanPSDAcrossTrialsElectrodesST{r,c,:};

                meanPSDAcrossTrialsLogElectrodesSessionsBL(r,c,:) = meanPSDAcrossTrialsLogElectrodesBL{r,c,:};
                meanPSDAcrossTrialsLogElectrodesSessionsST(r,c,:) = meanPSDAcrossTrialsLogElectrodesST{r,c,:};

                meanPSDAcrossLogTrialsElectrodesSessionsBL(r,c,:) = meanPSDAcrossLogTrialsElectrodesBL{r,c,:};
                meanPSDAcrossLogTrialsElectrodesSessionsST(r,c,:) = meanPSDAcrossLogTrialsElectrodesST{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hPSDFig(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hPSDFig(x,y));
            end
            if strcmpi(params.baselinestrategy,'conditionwise')
                plot(params.freqVals,squeeze(meanPSDAcrossTrialsLogElectrodesSessionsBL(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:),'linestyle','--');
            end
            set(gca,'nextplot','add');
            plot(params.freqVals,squeeze(meanPSDAcrossTrialsLogElectrodesSessionsST(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
            setAxesProperties(params.hPSDFig(x,y),params,x,y,rows);
            showTitleAndN(params.hPSDFig(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end
    set(params.hPSDFig,'ylim',ylimits);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean PSD across trials, electrodes and sessions','FitBoxToText','on');
    end    
    if strcmpi(params.baselinestrategy,'common')
        plot(params.freqVals,squeeze(conv2Log(mean(cell2mat([commonBaselinePSDAcrossTrials{:}]),2))),'color',[0.6 0.6 0.6],'linestyle','--');
    end

    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hPSDFig(r,c));
                delete(axisHandle);
            end
        end
    end 
end

%==========================================================================
% Difference PSDs: change of power from baseline to stimulus period
if ~isfield(params,'hdiffPSDFig') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hdiffPSDFig = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hdiffPSDFig')
    dPSDAcrossTrialsElectrodesSessions = zeros(numRows,numCols,numFreqPoints);
    dPSDAcrossTrialsLogElectrodesSessions = zeros(numRows,numCols,numFreqPoints);
    dPSDAcrossLogTrialsElectrodesSessions = zeros(numRows,numCols,numFreqPoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    dPSDAcrossTrialsElectrodesSessions(r,c,:) = squeeze((squeeze(cell2mat(dPSDAcrossTrialsElectrodes(r,c,:)))*weights)/sum(weights));
                    dPSDAcrossTrialsLogElectrodesSessions(r,c,:) = squeeze((squeeze(cell2mat(dPSDAcrossTrialsLogElectrodes(r,c,:)))*weights)/sum(weights));
                    dPSDAcrossLogTrialsElectrodesSessions(r,c,:) = squeeze((squeeze(cell2mat(dPSDAcrossLogTrialsElectrodes(r,c,:)))*weights)/sum(weights));
                else
                    dPSDAcrossTrialsElectrodesSessions(r,c,:) = squeeze(nanmean(squeeze(cell2mat(dPSDAcrossTrialsElectrodes(r,c,:))),2));
                    dPSDAcrossTrialsLogElectrodesSessions(r,c,:) = squeeze(nanmean(squeeze(cell2mat(dPSDAcrossTrialsLogElectrodes(r,c,:))),2));
                    dPSDAcrossLogTrialsElectrodesSessions(r,c,:) = squeeze(nanmean(squeeze(cell2mat(dPSDAcrossLogTrialsElectrodes(r,c,:))),2));
                end
            else
                dPSDAcrossTrialsElectrodesSessions(r,c,:) = dPSDAcrossTrialsElectrodes{r,c,:};
                dPSDAcrossTrialsLogElectrodesSessions(r,c,:) = dPSDAcrossTrialsLogElectrodes{r,c,:};
                dPSDAcrossLogTrialsElectrodesSessions(r,c,:) = dPSDAcrossLogTrialsElectrodes{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hdiffPSDFig(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hdiffPSDFig(x,y));
            end
            plot(params.freqVals,10*squeeze(dPSDAcrossTrialsLogElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:)); % in dB
            set(gca,'nextplot','add');
            line([0 params.fmax],[0 0],'color',[0.5 0.2 0.7]);
            setAxesProperties(params.hdiffPSDFig(x,y),params,x,y,rows);
            showTitleAndN(params.hdiffPSDFig(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end
    set(params.hdiffPSDFig,'ylim',ylimits);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean PSD across trials, electrodes and sessions','FitBoxToText','on');
    end
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hdiffPSDFig(r,c));
                delete(axisHandle);
            end
        end
    end 
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Difference PSDs: change of power from baseline to stimulus period
% With common baseline across conditions per electrode per session
if ~isfield(params,'hdiffPSDCommonBLFig') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hdiffPSDCommonBLFig = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hdiffPSDCommonBLFig')
dPSDAcrossTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols,numFreqPoints);
dPSDAcrossTrialsLogElectrodesSessionsCommonBL = zeros(numRows,numCols,numFreqPoints);
dPSDAcrossLogTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols,numFreqPoints);
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    dPSDAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze((squeeze(cell2mat(dPSDAcrossTrialsElectrodesCommonBL(r,c,:)))*weights)/sum(weights));
                    dPSDAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = squeeze((squeeze(cell2mat(dPSDAcrossTrialsLogElectrodesCommonBL(r,c,:)))*weights)/sum(weights));
                    dPSDAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze((squeeze(cell2mat(dPSDAcrossLogTrialsElectrodesCommonBL(r,c,:)))*weights)/sum(weights));
                else
                    dPSDAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(squeeze(cell2mat(dPSDAcrossTrialsElectrodesCommonBL(r,c,:))),2));
                    dPSDAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(squeeze(cell2mat(dPSDAcrossTrialsLogElectrodesCommonBL(r,c,:))),2));
                    dPSDAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(nanmean(squeeze(cell2mat(dPSDAcrossLogTrialsElectrodesCommonBL(r,c,:))),2));
                end
            else
                dPSDAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = dPSDAcrossTrialsElectrodesCommonBL{r,c,:};
                dPSDAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = dPSDAcrossTrialsLogElectrodesCommonBL{r,c,:};
                dPSDAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = dPSDAcrossLogTrialsElectrodesCommonBL{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hdiffPSDCommonBLFig(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hdiffPSDCommonBLFig(x,y));
            end
            plot(params.freqVals,10*squeeze(dPSDAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:)); % in dB
            set(gca,'nextplot','add');
            line([0 params.fmax],[0 0],'color',[0.5 0.2 0.7]);
            setAxesProperties(params.hdiffPSDCommonBLFig(x,y),params,x,y,rows);
            showTitleAndN(params.hdiffPSDCommonBLFig(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end
    set(params.hdiffPSDCommonBLFig,'ylim',ylimits);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean PSD across trials, electrodes and sessions','FitBoxToText','on');
    end
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hdiffPSDCommonBLFig(r,c));
                delete(axisHandle);
            end
        end
    end 
end

%==========================================================================
%--------------------------------------------------------------------------
% Plot marginals
%--------------------------------------------------------------------------
if isfield(params,'showmarginals')
    if params.showmarginals
        % Row-wise
        if ~isfield(params,'hPSDMarginalsFigRow') && newfiguresflag
            figure;
            plotsPos = [0.08 0.15 0.87 0.8];
            params.hPSDMarginalsFigRow = getPlotHandles(numRows,1,plotsPos);
        end
        
        if isfield(params,'hPSDMarginalsFigRow')
            for r = 1:numRows % one plot per row
                subplot(params.hPSDMarginalsFigRow(r,1));
                for c = 1:numCols
                    plot(params.freqVals,squeeze(meanPSDAcrossTrialsLogElectrodesSessionsBL(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:),'linestyle','--');
                    set(gca,'nextplot','add');
                    phandle{c} = plot(params.freqVals,squeeze(meanPSDAcrossTrialsLogElectrodesSessionsST(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
                    title(params.titleStringY{r,1});
                end

                if r==numRows
                    setAxesProperties(gca,params,1,1,1);
                end
                if r==ceil(numRows/2)
                    legend([phandle{:}],params.titleStringX);
                    title(params.titleStringY{r,1});
                end
            end
            if newfiguresflag
                annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
                annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean PSD across trials, electrodes and sessions','FitBoxToText','on');
            end
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % column-wise
        if ~isfield(params,'hPSDMarginalsFigCol') && newfiguresflag
            figure;
            plotsPos = [0.08 0.15 0.87 0.8];
            params.hPSDMarginalsFigCol = getPlotHandles(1,numCols,plotsPos);
        end
        if isfield(params,'hPSDMarginalsFigCol')
            clear phandle
            for c = 1:numCols % one plot per col
                subplot(params.hPSDMarginalsFigCol(1,c));
                for r = 1:numRows
                    plot(params.freqVals,squeeze(meanPSDAcrossTrialsLogElectrodesSessionsBL(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:),'linestyle','--');
                    set(gca,'nextplot','add');
                    phandle{r} = plot(params.freqVals,squeeze(meanPSDAcrossTrialsLogElectrodesSessionsST(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
                    title(params.titleStringX{1,c});
                end
                if c==1
                    setAxesProperties(gca,params,1,1,1);
                end
                if c==ceil(numCols/2)
                    legend([phandle{:}],params.titleStringY{:,c});
                end
            end
            if newfiguresflag
                annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
                annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean PSD across trials, electrodes and sessions','FitBoxToText','on');
            end
        end
    end
end
%==========================================================================
%--------------------------------------------------------------------------
% Plot Trends for the PSD measures (band power, peak frequency)
%--------------------------------------------------------------------------

if isfield(params,'showtrends') && isfield(params,'polar')
    dBandPowerAcrossTrialsElectrodesSessions = zeros(numRows,numCols);
    dBandPowerAcrossTrialsLogElectrodesSessions = zeros(numRows,numCols);
    dBandPowerAcrossLogTrialsElectrodesSessions = zeros(numRows,numCols);
    
    dBandPowerAcrossTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols);
    dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL = zeros(numRows,numCols);
    dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols);
    
    peakFreqdBandPowerAcrossTrialsElectrodesSessions = zeros(numRows,numCols);
    peakFreqdBandPowerAcrossTrialsLogElectrodesSessions = zeros(numRows,numCols);
    peakFreqdBandPowerAcrossLogTrialsElectrodesSessions = zeros(numRows,numCols);
    
    peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols);
    peakFreqdBandPowerAcrossTrialsLogElectrodesSessionsCommonBL = zeros(numRows,numCols);
    peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL = zeros(numRows,numCols);
    
    if params.showtrends || params.polar % these computations are done and plots are plotted when any one of 'trends' or 'polar' is selected (or both are selected)
        for r = 1:numRows
            for c = 1:numCols
                if numProtocols==1 % take average across protocols only when analyzing multiple protocols
                    dBandPowerAcrossTrialsElectrodesSessions(r,c,:) = dBandPowerAcrossTrialsElectrodes{r,c,:};
                    dBandPowerAcrossTrialsLogElectrodesSessions(r,c,:) = dBandPowerAcrossTrialsLogElectrodes{r,c,:};
                    dBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = dBandPowerAcrossLogTrialsElectrodes{r,c,:};

                    dBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = dBandPowerAcrossTrialsElectrodesCommonBL{r,c,:};
                    dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = dBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,:};
                    dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = dBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,:};

                    peakFreqdBandPowerAcrossTrialsElectrodesSessions(r,c,:) = peakFreqdBandPowerAcrossTrialsElectrodes{r,c,:};
                    peakFreqdBandPowerAcrossTrialsLogElectrodesSessions(r,c,:) = peakFreqdBandPowerAcrossTrialsLogElectrodes{r,c,:};
                    peakFreqdBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = peakFreqdBandPowerAcrossLogTrialsElectrodes{r,c,:};

                    peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = peakFreqdBandPowerAcrossTrialsElectrodesCommonBL{r,c,:};
                    peakFreqdBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = peakFreqdBandPowerAcrossTrialsLogElectrodesCommonBL{r,c,:};
                    peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL{r,c,:};
                else
                     if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                        clear weights
                        weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')'; % the total number of repeats across all electrodes for each session
                        dBandPowerAcrossTrialsElectrodesSessions(r,c,:) = squeeze(sum((squeeze(cell2mat(dBandPowerAcrossTrialsElectrodes(r,c,:))).*weights))./sum(weights));
                        dBandPowerAcrossTrialsLogElectrodesSessions(r,c,:) = squeeze(sum((squeeze(cell2mat(dBandPowerAcrossTrialsLogElectrodes(r,c,:))).*weights))./sum(weights));
                        dBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = squeeze(sum((squeeze(cell2mat(dBandPowerAcrossLogTrialsElectrodes(r,c,:))).*weights))./sum(weights));
                        % the value for a particular condition for each session
                        % is weighted by the respective number of repeats
                        % (across all electrodes considered) and the weighted
                        % mean is taken accordingly
                        dBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(sum((squeeze(cell2mat(dBandPowerAcrossTrialsElectrodesCommonBL(r,c,:))).*weights))./sum(weights));
                        dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = squeeze(sum((squeeze(cell2mat(dBandPowerAcrossTrialsLogElectrodesCommonBL(r,c,:))).*weights))./sum(weights));
                        dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(sum((squeeze(cell2mat(dBandPowerAcrossLogTrialsElectrodesCommonBL(r,c,:))).*weights))./sum(weights));

                        peakFreqdBandPowerAcrossTrialsElectrodesSessions(r,c,:) = squeeze(sum((squeeze(cell2mat(peakFreqdBandPowerAcrossTrialsElectrodes(r,c,:))).*weights))./sum(weights));
                        peakFreqdBandPowerAcrossTrialsLogElectrodesSessions(r,c,:) = squeeze(sum((squeeze(cell2mat(peakFreqdBandPowerAcrossTrialsLogElectrodes(r,c,:))).*weights))./sum(weights));
                        peakFreqdBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = squeeze(sum((squeeze(cell2mat(peakFreqdBandPowerAcrossLogTrialsElectrodes(r,c,:))).*weights))./sum(weights));

                        peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(sum((squeeze(cell2mat(peakFreqdBandPowerAcrossTrialsElectrodesCommonBL(r,c,:))).*weights))./sum(weights));
                        peakFreqdBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = squeeze(sum((squeeze(cell2mat(peakFreqdBandPowerAcrossTrialsLogElectrodesCommonBL(r,c,:))).*weights))./sum(weights));
                        peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = squeeze(sum((squeeze(cell2mat(peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL(r,c,:))).*weights))./sum(weights));
                     else
                        dBandPowerAcrossTrialsElectrodesSessions(r,c,:) = mean(squeeze(cell2mat(dBandPowerAcrossTrialsElectrodes(r,c,:))));
                        dBandPowerAcrossTrialsLogElectrodesSessions(r,c,:) = mean(squeeze(cell2mat(dBandPowerAcrossTrialsLogElectrodes(r,c,:))));
                        dBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = mean(squeeze(cell2mat(dBandPowerAcrossLogTrialsElectrodes(r,c,:))));

                        dBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = mean(squeeze(cell2mat(dBandPowerAcrossTrialsElectrodesCommonBL(r,c,:))));
                        dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = mean(squeeze(cell2mat(dBandPowerAcrossTrialsLogElectrodesCommonBL(r,c,:))));
                        dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = mean(squeeze(cell2mat(dBandPowerAcrossLogTrialsElectrodesCommonBL(r,c,:))));

                        peakFreqdBandPowerAcrossTrialsElectrodesSessions(r,c,:) = mean(squeeze(cell2mat(peakFreqdBandPowerAcrossTrialsElectrodes(r,c,:))));
                        peakFreqdBandPowerAcrossTrialsLogElectrodesSessions(r,c,:) = mean(squeeze(cell2mat(peakFreqdBandPowerAcrossTrialsLogElectrodes(r,c,:))));
                        peakFreqdBandPowerAcrossLogTrialsElectrodesSessions(r,c,:) = mean(squeeze(cell2mat(peakFreqdBandPowerAcrossLogTrialsElectrodes(r,c,:))));

                        peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL(r,c,:) = mean(squeeze(cell2mat(peakFreqdBandPowerAcrossTrialsElectrodesCommonBL(r,c,:))));
                        peakFreqdBandPowerAcrossTrialsLogElectrodesSessionsCommonBL(r,c,:) = mean(squeeze(cell2mat(peakFreqdBandPowerAcrossTrialsLogElectrodesCommonBL(r,c,:))));
                        peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL(r,c,:) = mean(squeeze(cell2mat(peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL(r,c,:))));
                     end
                end
            end
        end
    end
    
    if params.showtrends
        if ~isfield(params,'hPSDdBandPowerTrend') && newfiguresflag
            figure;
            plotsPos = [0.08 0.15 0.87 0.8];
            params.hPSDdBandPowerTrend = getPlotHandles(numRows,1,plotsPos);
        end
        
        if isfield(params,'hPSDdBandPowerTrend')
            % Plot change in power in the selected band
            subplot(params.hPSDdBandPowerTrend); 
            for r = 1:numRows
                for c = 1:numCols
                    plot(params.valsUnique{params.paramx}(c),squeeze(dBandPowerAcrossTrialsLogElectrodesSessions(r,c,:)),'o','color',params.colorsList((r-1)*numCols+c,:));
                    set(gca,'nextplot','add');
                end
                plot(params.valsUnique{params.paramx},squeeze(dBandPowerAcrossTrialsLogElectrodesSessions(r,:,:)),'color',[0.8 0.8 0.8]/r);
                linename = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
                text(0.1,0.95-0.1*r,linename,'fontsize',6,'color',[0.8 0.8 0.8]/r,'unit','normalized');
            end
            title(['Change in power in the frequency band [' num2str(params.fBandLow) ' ' num2str(params.fBandHigh) '] Hz']);
            ylabel('log(P_s/P_b)');
            line([0 params.valsUnique{params.paramx}(numCols)],[0 0],'color',[0.2 0.7 0.4]);
            setAxesBasicProperties(gca);
        end
        
        if ~isfield(params,'hPSDPeakFreqTrend') && newfiguresflag
            figure;
            plotsPos = [0.08 0.15 0.87 0.8];
            params.hPSDPeakFreqTrend = getPlotHandles(numRows,1,plotsPos);
        end
        
        if isfield(params,'hPSDPeakFreqTrend')
            % Plot peak frequency in the selected band
            subplot(params.hPSDPeakFreqTrend);
            for r = 1:numRows
                for c = 1:numCols
                    plot(params.valsUnique{params.paramx}(c),squeeze(peakFreqdBandPowerAcrossTrialsLogElectrodesSessions(r,c,:)),'o','color',params.colorsList((r-1)*numCols+c,:));
                    set(gca,'nextplot','add');
                end
                plot(params.valsUnique{params.paramx},squeeze(peakFreqdBandPowerAcrossTrialsLogElectrodesSessions(r,:,:)),'color',[0.8 0.8 0.8]/r);
                linename = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
                text(0.1,0.95-0.1*r,linename,'fontsize',6,'color',[0.8 0.8 0.8]/r,'unit','normalized');
            end
            title(['Peak frequency in the frequency band [' num2str(params.fBandLow) ' ' num2str(params.fBandHigh) '] Hz']);
            ylabel('Frequency (Hz)');
            line([0 params.valsUnique{params.paramx}(numCols)],[0 0],'color',[0.2 0.7 0.4]);
            setAxesBasicProperties(gca);
            set(gca,'ylim',[params.fBandLow params.fBandHigh]);
        end
    end
end

%--------------------------------------------------------------------------
% Make polar plots
%--------------------------------------------------------------------------
if isfield(params,'polar')
    if params.polar
        % check which param corresponds to orientation and assign numOris
        % accrodingly
        if params.paramx==5
            numOris = numCols;
            numCases = numRows;
            polarcase = 1;
        elseif params.paramy==5
            numOris = numRows;
            numCases = numCols;
            polarcase = 2;
        else
            numOris = numCols;
            numCases = numRows;
            polarcase = 1;
        end

        % Ask if there are any non hue stimuli (eg. grating) for which only the
        % corresponding data point will be shown at 0 deg hue
        extraDataPoints = str2double(inputdlg('Are there any non-hue stimuli or pooled conditions? How many?'));

        if polarcase==1 % cols are represented as hues
            for r=1:numCases
                if ~isfield(params,'hPSDBandPowerPolar') && newfiguresflag
                    handlesDontExist = 1;
                else
                    handlesDontExist = 0;
                end
                if handlesDontExist
                    figure;
                    params.hPSDBandPowerPolar = subplot(2,2,1);
                end
                filldata = squeeze(dBandPowerAcrossTrialsLogElectrodesSessions(r,:,:));
                numTheta = numOris-extraDataPoints;
                hValsList = params.valsUnique{5}(1:numTheta); % ori = hue
                if extraDataPoints
                    extraPts = filldata(end-extraDataPoints:end);
                    makePolarPlotHSV(params.hPSDBandPowerPolar,filldata,numTheta,hValsList,'linear','extra',num2cell(extraPts));
                else
                    makePolarPlotHSV(params.hPSDBandPowerPolar,filldata,numTheta,hValsList,'linear');
                end
                title('PSD band power');
                setAxesBasicProperties(params.hPSDBandPowerPolar);
                set(params.hPSDBandPowerPolar,'yticklabel',[]);
                set(params.hPSDBandPowerPolar,'xticklabel',[]);

                if handlesDontExist
                    polarfigname = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
                    set(gcf,'numbertitle','off','name',polarfigname);
                end
            end

        else % rows are represented as hues
            for c=1:numCases
                if ~isfield(params,'hPSDBandPowerPolar') && newfiguresflag
                    handlesDontExist = 1;
                else
                    handlesDontExist = 0;
                end
                if handlesDontExist
                    figure;
                    params.hPSDBandPowerPolar = subplot(2,2,1);
                end

                polarfigname = [params.paramsString{params.paramx} ':' num2str(params.valsUnique{params.paramx}(c))];
                filldata = (squeeze(meanstimPeriodMin(:,c,:)))';
                numTheta = numOris-extraDataPoints;
                hValsList = params.valsUnique{5}(1:numTheta); % ori = hue
                makePolarPlotHSV(params.hPSDBandPowerPolar,filldata,numTheta,hValsList,'abs','extra',extraDataPoints);
                title('PSD Min');
                setAxesBasicProperties(params.hPSDBandPowerPolar);
                set(params.hPSDBandPowerPolar,'yticklabel',[]);
                set(params.hPSDBandPowerPolar,'xticklabel',[]);

                set(gcf,'numbertitle','off','name',polarfigname);
            end
        end
    end
end

%==========================================================================
%--------------------------------------------------------------------------
% Make topoplots
%--------------------------------------------------------------------------
if isfield(params,'topo')
    dBandPowerPerElectrodeAcrossTrialsSessions = zeros(numRows,numElecs,numCols);
    dBandPowerPerElectrodeAcrossLogTrialsSessions = zeros(numRows,numElecs,numCols);
    peakFreqdBandPowerPerElectrodeAcrossTrialsSessions = zeros(numRows,numElecs,numCols);
    peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessions = zeros(numRows,numElecs,numCols);
    if params.topo
        for r=1:numRows
            if numProtocols>1 % take average across protocols only when analyzing multiple protocols
                if params.weightedmean && strcmpi(params.badtrialsoption,'common')
                    clear weights
                    for c=1:numCols
                        weights = cell2mat(squeeze(params.numRepeats(r,c,:,:))); % contains the repeats per electrode per session
                        % weights: numCols x numProtocols
                        dBandPowerPerElectrodeAcrossTrialsSessions(r,:,c) = sum((weights.*(reshape(cell2mat([dBandPowerAcrossTrials{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                        peakFreqdBandPowerPerElectrodeAcrossTrialsSessions(r,:,c) = sum((weights.*(reshape(cell2mat([peakFreqdBandPowerAcrossTrials{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);

                        dBandPowerPerElectrodeAcrossLogTrialsSessions(r,:,c) = sum((weights.*(reshape(cell2mat([dBandPowerAcrossLogTrials{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                        peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessions(r,:,c) = sum((weights.*(reshape(cell2mat([peakFreqdBandPowerAcrossLogTrials{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);                    
                    end
                elseif params.weightedmean && ~strcmpi(params.badtrialsoption,'common') % individual bad trials
                    clear weights
                    for c=1:numCols
                        weights = cell2mat(squeeze(params.numRepeats(r,c,:,:)));
                        % weight: numProtocols x numElecs 
                        % We use sum here because it is not a matrix multiplication
                        dBandPowerPerElectrodeAcrossTrialsSessions(r,:,c) = sum((weights.*(reshape(cell2mat([dBandPowerAcrossTrials{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                        peakFreqdBandPowerPerElectrodeAcrossTrialsSessions(r,:,c) = sum((weights.*(reshape(cell2mat([peakFreqdBandPowerAcrossTrials{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);

                        dBandPowerPerElectrodeAcrossLogTrialsSessions(r,:,c) = sum((weights.*(reshape(cell2mat([dBandPowerAcrossLogTrials{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                        peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessions(r,:,c) = sum((weights.*(reshape(cell2mat([peakFreqdBandPowerAcrossLogTrials{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                    end
                else
                    dBandPowerPerElectrodeAcrossTrialsSessions(r,:,c) = nanmean(reshape(cell2mat([dBandPowerAcrossTrials{r,c,:}]),numElecs,numProtocols),2);
                    peakFreqdBandPowerPerElectrodeAcrossTrialsSessions(r,:,c) = nanmean(reshape(cell2mat([peakFreqdBandPowerAcrossTrials{r,c,:}]),numElecs,numProtocols),2);

                    dBandPowerPerElectrodeAcrossLogTrialsSessions(r,:,c) = nanmean(reshape(cell2mat([dBandPowerAcrossLogTrials{r,c,:}]),numElecs,numProtocols),2);
                    peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessions(r,:,c) = nanmean(reshape(cell2mat([peakFreqdBandPowerAcrossLogTrials{r,c,:}]),numElecs,numProtocols),2);
                end

            else
                for c=1:numCols
                    dBandPowerPerElectrodeAcrossTrialsSessions(r,:,c) = cell2mat(dBandPowerAcrossTrials{r,c,:});
                    peakFreqdBandPowerPerElectrodeAcrossTrialsSessions(r,:,c) = cell2mat(peakFreqdBandPowerAcrossTrials{r,c,:});

                    dBandPowerPerElectrodeAcrossLogTrialsSessions(r,:,c) = cell2mat(dBandPowerAcrossLogTrials{r,c,:});
                    peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessions(r,:,c) = cell2mat(peakFreqdBandPowerAcrossLogTrials{r,c,:});
                end
            end
        end

        % decide the number of topoplots and assign data accordingly
        clear dataAllConditions

        dataAllConditions{1} = dBandPowerPerElectrodeAcrossTrialsSessions;
        dataAllConditions{2} = peakFreqdBandPowerPerElectrodeAcrossTrialsSessions;
        figname{1} = 'topoplot_bandPowerPSD';
        figname{2} = 'topoplot_peakFreqPSD';

        % Draw the topoplots
        handlesDontExist = 0;
        for k = 1:length(dataAllConditions)  
            if ~isfield(params,'hTopoPSDFig') && newfiguresflag
                handlesDontExist = 1;
            end
            if handlesDontExist
                figure;
                plotsPos = [0.08 0.15 0.87 0.8];
                if numCols>params.maxPlotsAlongX
                    cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
                    subplotsRearrangement = 1;
                else
                    cols = numCols; rows = numRows; subplotsRearrangement = 0;
                end
                params.hTopoPSDFig{k} = getPlotHandles(rows,cols,plotsPos);
            end

            if ~iscell(params.hTopoPSDFig)
                tempArray = params.hTopoPSDFig;
                params = rmfield(params,'hTopoPSDFig');
                params.hTopoPSDFig{1} = tempArray;
            end

            colormap(jet);

            for r = 1:numRows
                for c = 1:numCols
                    if ~subplotsRearrangement
                        x = r; y = c;
                        subplot(params.hTopoPSDFig{k}(x,y));
                    else
                        x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                        y = mod(c,cols);
                        if y==0; y = cols; end
                        subplot(params.hTopoPSDFig{k}(x,y));
                    end

                    datatoplot = zeros(1,length(params.chanlocs));
                    datatoplot(:,params.electrodes) = squeeze(dataAllConditions{k}(r,:,c));
                    topoplot(datatoplot,params.chanlocs,'electrodes','off','drawaxis','off',...
                        'plotchans',params.electrodes,'hcolor',[0.25 0.25 0.25]);
                    whitebg([0 0 0]);
                    showTitleAndN(params.hTopoPSDFig{k}(x,y),params,r,c,0,-0.25,0.94); % show title and n outside the top left part of topoplot
                    drawnow;
                end
            end
            set(params.hTopoPSDFig{k},'clim',[min(min(min(dataAllConditions{k}))) max(max(max(dataAllConditions{k})))]);
            if newfiguresflag
                annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
                annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean across trials and sessions','FitBoxToText','on');
                colorbar('position',[plotsPos(1)+plotsPos(3)+0.005 plotsPos(2) 0.015 plotsPos(4)]);
                set(gcf,'numbertitle','off','name',figname{k});
            end
            % delete unused axes
            for r=1:rows
                for c = 1:cols
                    if cols*(r-1)+c>numRows*numCols
                        axisHandle = findall(params.hTopoPSDFig{k}(r,c));
                        delete(axisHandle);
                    end
                end
            end 
        end

        %______________________________________________________________________
        % with common baseline
        dBandPowerPerElectrodeAcrossTrialsSessionsCommonBL = zeros(numRows,numElecs,numCols);
        dBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL = zeros(numRows,numElecs,numCols);
        peakFreqdBandPowerPerElectrodeAcrossTrialsSessionsCommonBL = zeros(numRows,numElecs,numCols);
        peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL = zeros(numRows,numElecs,numCols);
        for r=1:numRows
            if numProtocols>1 % take average across protocols only when analyzing multiple protocols
                if params.weightedmean && strcmpi(params.badtrialsoption,'common')
                    clear weights
                    for c=1:numCols
                        weights = cell2mat(squeeze(params.numRepeats(r,c,:,:))); % contains the repeats per electrode per session
                        % weights: numCols x numProtocols
                        dBandPowerPerElectrodeAcrossTrialsSessionsCommonBL(r,:,c) = sum((weights.*(reshape(cell2mat([dBandPowerAcrossTrialsCommonBL{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                        peakFreqdBandPowerPerElectrodeAcrossTrialsSessionsCommonBL(r,:,c) = sum((weights.*(reshape(cell2mat([peakFreqdBandPowerAcrossTrialsCommonBL{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);

                        dBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL(r,:,c) = sum((weights.*(reshape(cell2mat([dBandPowerAcrossLogTrialsCommonBL{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                        peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL(r,:,c) = sum((weights.*(reshape(cell2mat([peakFreqdBandPowerAcrossLogTrialsCommonBL{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                    end
                elseif params.weightedmean && ~strcmpi(params.badtrialsoption,'common') % individual bad trials
                    clear weights
                    for c=1:numCols
                        weights = cell2mat(squeeze(params.numRepeats(r,c,:,:)));
                        % weight: numProtocols x numElecs 
                        % We use sum here because it is not a matrix multiplication
                        dBandPowerPerElectrodeAcrossTrialsSessionsCommonBL(r,:,c) = sum((weights.*(reshape(cell2mat([dBandPowerAcrossTrialsCommonBL{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                        peakFreqdBandPowerPerElectrodeAcrossTrialsSessionsCommonBL(r,:,c) = sum((weights.*(reshape(cell2mat([peakFreqdBandPowerAcrossTrialsCommonBL{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);

                        dBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL(r,:,c) = sum((weights.*(reshape(cell2mat([dBandPowerAcrossLogTrialsCommonBL{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                        peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL(r,:,c) = sum((weights.*(reshape(cell2mat([peakFreqdBandPowerAcrossLogTrialsCommonBL{r,c,:}]),numElecs,numProtocols))'),1)./sum(weights,1);
                    end
                else
                    dBandPowerPerElectrodeAcrossTrialsSessionsCommonBL(r,:,c) = nanmean(reshape(cell2mat([dBandPowerAcrossTrialsCommonBL{r,c,:}]),numElecs,numProtocols),2);
                    peakFreqdBandPowerPerElectrodeAcrossTrialsSessionsCommonBL(r,:,c) = nanmean(reshape(cell2mat([peakFreqdBandPowerAcrossTrialsCommonBL{r,c,:}]),numElecs,numProtocols),2);

                    dBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL(r,:,c) = nanmean(reshape(cell2mat([dBandPowerAcrossLogTrialsCommonBL{r,c,:}]),numElecs,numProtocols),2);
                    peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL(r,:,c) = nanmean(reshape(cell2mat([peakFreqdBandPowerAcrossLogTrialsCommonBL{r,c,:}]),numElecs,numProtocols),2);
                end

            else
                for c=1:numCols
                    dBandPowerPerElectrodeAcrossTrialsSessionsCommonBL(r,:,c) = cell2mat(dBandPowerAcrossTrialsCommonBL{r,c,:});
                    peakFreqdBandPowerPerElectrodeAcrossTrialsSessionsCommonBL(r,:,c) = cell2mat(peakFreqdBandPowerAcrossTrialsCommonBL{r,c,:});

                    dBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL(r,:,c) = cell2mat(dBandPowerAcrossLogTrialsCommonBL{r,c,:});
                    peakFreqdBandPowerPerElectrodeAcrossLogTrialsSessionsCommonBL(r,:,c) = cell2mat(peakFreqdBandPowerAcrossLogTrialsCommonBL{r,c,:});
                end
            end
        end

        % decide the number of topoplots and assign data accordingly
        clear dataAllConditions

        dataAllConditions{1} = dBandPowerPerElectrodeAcrossTrialsSessionsCommonBL;
        dataAllConditions{2} = peakFreqdBandPowerPerElectrodeAcrossTrialsSessionsCommonBL;
        figname{1} = 'topoplot_bandPowerPSD';
        figname{2} = 'topoplot_peakFreqPSD';

        % Draw the topoplots
        handlesDontExist = 0;
        for k = 1:length(dataAllConditions)  
            if ~isfield(params,'hTopoPSDCommonBLFig') && newfiguresflag
                handlesDontExist = 1;
            end
            if handlesDontExist
                figure;
                plotsPos = [0.08 0.15 0.87 0.8];
                if numCols>params.maxPlotsAlongX
                    cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
                    subplotsRearrangement = 1;
                else
                    cols = numCols; rows = numRows; subplotsRearrangement = 0;
                end
                params.hTopoPSDCommonBLFig{k} = getPlotHandles(rows,cols,plotsPos);
            end

            if ~iscell(params.hTopoPSDCommonBLFig)
                tempArray = params.hTopoPSDCommonBLFig;
                params = rmfield(params,'hTopoPSDFigCommonBL');
                params.hTopoPSDCommonBLFig{1} = tempArray;
            end

            colormap(jet);

            for r = 1:numRows
                for c = 1:numCols
                    if ~subplotsRearrangement
                        x = r; y = c;
                        subplot(params.hTopoPSDCommonBLFig{k}(x,y));
                    else
                        x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                        y = mod(c,cols);
                        if y==0; y = cols; end
                        subplot(params.hTopoPSDCommonBLFig{k}(x,y));
                    end

                    datatoplot = zeros(1,length(params.chanlocs));
                    datatoplot(:,params.electrodes) = squeeze(dataAllConditions{k}(r,:,c));
                    topoplot(datatoplot,params.chanlocs,'electrodes','off','drawaxis','off',...
                        'plotchans',params.electrodes,'hcolor',[0.25 0.25 0.25]);
                    whitebg([0 0 0]);
                    showTitleAndN(params.hTopoPSDCommonBLFig{k}(x,y),params,r,c,0,-0.25,0.94); % show title and n outside the top left part of topoplot
                    drawnow;
                end
            end
            set(params.hTopoPSDCommonBLFig{k},'clim',[min(min(min(dataAllConditions{k}))) max(max(max(dataAllConditions{k})))]);
            if newfiguresflag
                annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
                annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean across trials and sessions','FitBoxToText','on');
                colorbar('position',[plotsPos(1)+plotsPos(3)+0.005 plotsPos(2) 0.015 plotsPos(4)]);
                set(gcf,'numbertitle','off','name',figname{k});
            end
            % delete unused axes
            for r=1:rows
                for c = 1:cols
                    if cols*(r-1)+c>numRows*numCols
                        axisHandle = findall(params.hTopoPSDCommonBLFig{k}(r,c));
                        delete(axisHandle);
                    end
                end
            end 
        end    
    end
end
%--------------------------------------------------------------------------
end
%==========================================================================
%==========================================================================
% return the specified psd data
%--------------------------------------------------------------------------
if returnData
    if strcmpi(params.baselinestrategy,'conditionwise')
        switch params.psdstrategy
            case 'raw' % no logarithmic operation
                psddata.psdBL = PSDAcrossTrialsBL;
                psddata.psdST = PSDAcrossTrialsST;
                psddata.psdMeanBL = meanPSDAcrossTrialsElectrodesBL;
                psddata.psdMeanST = meanPSDAcrossTrialsElectrodesST;
                psddata.diffPSD = dPSDAcrossTrials;
                psddata.diffPSDMean = dPSDAcrossTrialsElectrodes;
                psddata.dBandPower = dBandPowerAcrossTrials;
                psddata.dBandPowerMean = dBandPowerAcrossTrialsElectrodes;
                psddata.peakFreq = peakFreqdBandPowerAcrossTrials;
                psddata.peakFreqMean = peakFreqdBandPowerAcrossTrialsElectrodes;
                if drawplots && newfiguresflag
                    psddata.psdMeanSessionsBL = meanPSDAcrossTrialsElectrodesSessionsBL;
                    psddata.psdMeanSessionsST = meanPSDAcrossTrialsElectrodesSessionsST;
                    psddata.diffPSDMeanSessions = dPSDAcrossTrialsElectrodesSessions;
                    psddata.dBandPowerMeanSessions = dBandPowerAcrossTrialsElectrodesSessions;
                    psddata.peakFreqMeanSessions = peakFreqdBandPowerAcrossTrialsElectrodesSessions;
                end
                if isfield(params,'extraFreqBands')
                    psddata.dBandPowerEB = dBandPowerAcrossTrialsEB;
                    psddata.dBandPowerMeanEB = dBandPowerAcrossTrialsElectrodesEB;
                    psddata.peakFreqEB = peakFreqdBandPowerAcrossTrialsEB;
                    psddata.peakFreqMeanEB = peakFreqdBandPowerAcrossTrialsElectrodesEB;
                end

            case 'logTrial' % take log PSD of every trial
                psddata.psdBL = PSDAcrossLogTrialsBL;
                psddata.psdST = PSDAcrossLogTrialsST;
                psddata.psdMeanBL = meanPSDAcrossLogTrialsElectrodesBL;
                psddata.psdMeanST = meanPSDAcrossLogTrialsElectrodesST;
                psddata.diffPSD = dPSDAcrossLogTrials;
                psddata.diffPSDMean = dPSDAcrossLogTrialsElectrodes;
                psddata.dBandPower = dBandPowerAcrossLogTrials;
                psddata.dBandPowerMean = dBandPowerAcrossLogTrialsElectrodes;
                psddata.peakFreq = peakFreqdBandPowerAcrossLogTrials;
                psddata.peakFreqMean = peakFreqdBandPowerAcrossLogTrialsElectrodes;
                if drawplots && newfiguresflag
                    psddata.psdMeanSessionsBL = meanPSDAcrossLogTrialsElectrodesSessionsBL;
                    psddata.psdMeanSessionsST = meanPSDAcrossLogTrialsElectrodesSessionsST;
                    psddata.diffPSDMeanSessions = dPSDAcrossLogTrialsElectrodesSessions;
                    psddata.dBandPowerMeanSessions = dBandPowerAcrossLogTrialsElectrodesSessions;
                    psddata.peakFreqMeanSessions = peakFreqdBandPowerAcrossLogTrialsElectrodesSessions;
                end
                if isfield(params,'extraFreqBands')
                    psddata.dBandPowerEB = dBandPowerAcrossLogTrialsEB;
                    psddata.dBandPowerMeanEB = dBandPowerAcrossLogTrialsElectrodesEB;
                    psddata.peakFreqEB = peakFreqdBandPowerAcrossLogTrialsEB;
                    psddata.peakFreqMeanEB = peakFreqdBandPowerAcrossLogTrialsElectrodesEB;
                end
            otherwise % default - take log at the level of every electrode
                psddata.psdBL = PSDAcrossTrialsBL;
                psddata.psdST = PSDAcrossTrialsST;
                psddata.psdMeanBL = meanPSDAcrossTrialsLogElectrodesBL;
                psddata.psdMeanST = meanPSDAcrossTrialsLogElectrodesST;
                psddata.diffPSD = dPSDAcrossTrials;
                psddata.diffPSDMean = dPSDAcrossTrialsLogElectrodes;
                psddata.dBandPower = dBandPowerAcrossTrials;
                psddata.dBandPowerMean = dBandPowerAcrossTrialsLogElectrodes;
                psddata.peakFreq = peakFreqdBandPowerAcrossTrials;
                psddata.peakFreqMean = peakFreqdBandPowerAcrossTrialsElectrodes;
                if drawplots && newfiguresflag
                    psddata.psdMeanSessionsBL = meanPSDAcrossTrialsLogElectrodesSessionsBL;
                    psddata.psdMeanSessionsST = meanPSDAcrossTrialsLogElectrodesSessionsST;
                    psddata.diffPSDMeanSessions = dPSDAcrossTrialsLogElectrodesSessions;
                    psddata.dBandPowerMeanSessions = dBandPowerAcrossTrialsLogElectrodesSessions;
                    psddata.peakFreqMeanSessions = peakFreqdBandPowerAcrossTrialsLogElectrodesSessions;
                end
                if isfield(params,'extraFreqBands')
                    psddata.dBandPowerEB = dBandPowerAcrossTrialsEB;
                    psddata.dBandPowerMeanEB = dBandPowerAcrossTrialsLogElectrodesEB;
                    psddata.peakFreqEB = peakFreqdBandPowerAcrossTrialsEB;
                    psddata.peakFreqMeanEB = peakFreqdBandPowerAcrossTrialsLogElectrodesEB;
                end
        end
    elseif strcmpi(params.baselinestrategy,'common')
        switch params.psdstrategy
            case 'raw' % no logarithmic operation
                psddata.psdBL = PSDAcrossTrialsBL;
                psddata.psdST = PSDAcrossTrialsST;
                psddata.psdMeanBL = meanPSDAcrossTrialsElectrodesBL;
                psddata.psdMeanST = meanPSDAcrossTrialsElectrodesST;
                psddata.diffPSD = dPSDAcrossTrialsCommonBL;
                psddata.diffPSDMean = dPSDAcrossTrialsElectrodesCommonBL;
                psddata.dBandPower = dBandPowerAcrossTrialsCommonBL;
                psddata.dBandPowerMean = dBandPowerAcrossTrialsElectrodesCommonBL;
                psddata.peakFreq = peakFreqdBandPowerAcrossTrialsCommonBL;
                psddata.peakFreqMean = peakFreqdBandPowerAcrossTrialsElectrodesCommonBL;
                if drawplots && newfiguresflag
                    psddata.psdMeanSessionsBL = meanPSDAcrossTrialsElectrodesSessionsBL;
                    psddata.psdMeanSessionsST = meanPSDAcrossTrialsElectrodesSessionsST;
                    psddata.diffPSDMeanSessions = dPSDAcrossTrialsElectrodesSessionsCommonBL;
                    psddata.dBandPowerMeanSessions = dBandPowerAcrossTrialsElectrodesSessionsCommonBL;
                    psddata.peakFreqMeanSessions = peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL;
                end
                if isfield(params,'extraFreqBands')
                    psddata.dBandPowerEB = dBandPowerAcrossTrialsCommonBLEB;
                    psddata.dBandPowerMeanEB = dBandPowerAcrossTrialsElectrodesCommonBLEB;
                    psddata.peakFreqEB = peakFreqdBandPowerAcrossTrialsCommonBLEB;
                    psddata.peakFreqMeanEB = peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB;
                end
                
                psddata.commonBaselinePSDTrials = commonBaselinePSDAcrossTrials;
                psddata.commonBaselinePSDLogTrials = commonBaselinePSDAcrossLogTrials;

            case 'logTrial' % take log PSD of every trial
                psddata.psdBL = PSDAcrossTrialsBL;
                psddata.psdST = PSDAcrossTrialsST;
                psddata.psdMeanBL = meanPSDAcrossLogTrialsElectrodesBL;
                psddata.psdMeanST = meanPSDAcrossLogTrialsElectrodesST;
                psddata.diffPSD = dPSDAcrossLogTrialsCommonBL;
                psddata.diffPSDMean = dPSDAcrossLogTrialsElectrodesCommonBL;
                psddata.dBandPower = dBandPowerAcrossLogTrialsCommonBL;
                psddata.dBandPowerMean = dBandPowerAcrossLogTrialsElectrodesCommonBL;
                psddata.peakFreq = peakFreqdBandPowerAcrossLogTrialsCommonBL;
                psddata.peakFreqMean = peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBL;
                if drawplots && newfiguresflag
                    psddata.psdMeanSessionsBL = meanPSDAcrossLogTrialsElectrodesSessionsBL;
                    psddata.psdMeanSessionsST = meanPSDAcrossLogTrialsElectrodesSessionsST;
                    psddata.diffPSDMeanSessions = dPSDAcrossLogTrialsElectrodesSessionsCommonBL;
                    psddata.dBandPowerMeanSessions = dBandPowerAcrossLogTrialsElectrodesSessionsCommonBL;
                    psddata.peakFreqMeanSessions = peakFreqdBandPowerAcrossLogTrialsElectrodesSessionsCommonBL;
                end
                if isfield(params,'extraFreqBands')
                    psddata.dBandPowerEB = dBandPowerAcrossLogTrialsCommonBLEB;
                    psddata.dBandPowerMeanEB = dBandPowerAcrossLogTrialsElectrodesCommonBLEB;
                    psddata.peakFreqEB = peakFreqdBandPowerAcrossLogTrialsCommonBLEB;
                    psddata.peakFreqMeanEB = peakFreqdBandPowerAcrossLogTrialsElectrodesCommonBLEB;
                end
                psddata.commonBaselinePSDTrials = commonBaselinePSDAcrossTrials;
                psddata.commonBaselinePSDLogTrials = commonBaselinePSDAcrossLogTrials;
                
            otherwise % default - take log at the level of every electrode
                psddata.psdBL = PSDAcrossTrialsBL;
                psddata.psdST = PSDAcrossTrialsST;
                psddata.psdMeanBL = meanPSDAcrossTrialsLogElectrodesBL;
                psddata.psdMeanST = meanPSDAcrossTrialsLogElectrodesST;
                psddata.diffPSD = dPSDAcrossTrialsCommonBL;
                psddata.diffPSDMean = dPSDAcrossTrialsLogElectrodesCommonBL;
                psddata.dBandPower = dBandPowerAcrossTrialsCommonBL;
                psddata.dBandPowerMean = dBandPowerAcrossTrialsLogElectrodesCommonBL;
                psddata.peakFreq = peakFreqdBandPowerAcrossTrialsCommonBL;
                psddata.peakFreqMean = peakFreqdBandPowerAcrossTrialsElectrodesCommonBL;
                if drawplots && newfiguresflag
                    psddata.psdMeanSessionsBL = meanPSDAcrossTrialsLogElectrodesSessionsBL;
                    psddata.psdMeanSessionsST = meanPSDAcrossTrialsLogElectrodesSessionsST;
                    psddata.diffPSDMeanSessions = dPSDAcrossTrialsLogElectrodesSessionsCommonBL;
                    psddata.dBandPowerMeanSessions = dBandPowerAcrossTrialsLogElectrodesSessionsCommonBL;
                    psddata.peakFreqMeanSessions = peakFreqdBandPowerAcrossTrialsElectrodesSessionsCommonBL;
                end
                if isfield(params,'extraFreqBands')
                    psddata.dBandPowerEB = dBandPowerAcrossTrialsCommonBLEB;
                    psddata.dBandPowerMeanEB = dBandPowerAcrossTrialsLogElectrodesCommonBLEB;
                    psddata.peakFreqEB = peakFreqdBandPowerAcrossTrialsCommonBLEB;
                    psddata.peakFreqMeanEB = peakFreqdBandPowerAcrossTrialsElectrodesCommonBLEB;
                end
                psddata.commonBaselinePSDTrials = commonBaselinePSDAcrossTrials;
                psddata.commonBaselinePSDLogTrials = commonBaselinePSDAcrossLogTrials;
        end
    end
    psddata.freqVals = params.freqVals;
    psddata.params = params;
else
    psddata = [];
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
    set(hAxes,'xlim',[params.fmin params.fmax]);
    set(hAxes,'tickLength',[0.03 0.01]);
    set(hAxes,'fontsize',20,'fontweight','bold');
    % set(hAxes,'xtick',[0 0.8],'tickdir','out');
    set(hAxes,'tickdir','out');
    set(hAxes,'box','off');

    if col==1 && row==rows
        ylabel(hAxes,'logPSD (uV^2)');
        xlabel(hAxes,'Frequency (Hz)');
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
set(hAxes,'xlim',[params.fmin params.fmax]);
set(hAxes,'tickLength',[0.03 0.01]);
set(hAxes,'fontsize',10,'fontweight','bold'); % use fontsize ~20 for presentation figures
set(hAxes,'tickdir','out');
set(hAxes,'box','off');  
end

function setAxesPropertiesFigure(hFig)
allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'xlim',[-0.5 1]);
set(allAxesInFigure,'tickLength',[0.03 0.01]);
set(allAxesInFigure,'fontsize',10,'fontweight','bold'); % use fontsize ~20 for presentation figures
set(allAxesInFigure,'xtick',[0 0.8],'tickdir','out');
set(allAxesInFigure,'box','off');
end