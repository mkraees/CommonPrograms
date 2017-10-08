% plot mean Firing rate (FR) across electrodes, sessions, conditions
% Vinay Shirhatti, 26 October 2016
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
% trialwise spikes for that electrode. For eg. if there are 16 trials for
% any particular condition for electrode number 3 in the list then:
% data{r,c,i}{3} : (16 x numTimeSamples) double array
%
% params : structure of variables/parameters that decide the computation
% drawplots : set this to 1 to draw the plots, otheriwse to 0 to only do
%             the calculations and return the FR data
% newfiguresflag: set this flag to 1 to plot everything. If the relevant 
%                 handles are passed then the plots are drawn in those
%                 subplots. Otherwise new figures are created and the plots
%                 are drawn in them.
% returnData : set this flag to 1 to return the relevant data, else frdata
%              is empty
%
% OUTPUTS
% params : modified params
% frdata : the calculated FR data
% *************************************************************************

function [params,frdata] = computeAndPlotMeanFRMeasures(data,params,drawplots,newfiguresflag,returnData,varargin)

if ~exist('drawplots','var'); drawplots=1; end
if ~exist('newfiguresflag','var'); newfiguresflag=0; end
if ~exist('returnData','var'); returnData=1; end

if sum(strcmpi('psth',varargin))
    psth = varargin{find(strcmpi(varargin,'psth'))+1};
end

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

if ~isfield(params,'FRsmooothSigmaMS')
    params.FRsmooothSigmaMS = []; % this variable stores the maximum number of plots to be drawn along the horizontal axis
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

%--------------------------------------------------------------------------
% Compute and plot FRs
%--------------------------------------------------------------------------
% Initialize variables
if ~exist('psth','var')
    psth = cell(numRows,numCols,numProtocols);
    computepsth = 1;
else
    computepsth = 0;
end
FRPerTrialBL = cell(numRows,numCols,numProtocols);
FRPerTrialST = cell(numRows,numCols,numProtocols);
FRAcrossTrialsBL = cell(numRows,numCols,numProtocols);
FRAcrossTrialsST = cell(numRows,numCols,numProtocols);
semFRAcrossTrialsBL = cell(numRows,numCols,numProtocols);
semFRAcrossTrialsST = cell(numRows,numCols,numProtocols);
meanFRAcrossTrialsElectrodesBL = cell(numRows,numCols,numProtocols);
meanFRAcrossTrialsElectrodesST = cell(numRows,numCols,numProtocols);
semFRAcrossTrialsElectrodesBL = cell(numRows,numCols,numProtocols);
semFRAcrossTrialsElectrodesST = cell(numRows,numCols,numProtocols);

psthNormalized = cell(numRows,numCols,numProtocols);

psthMeanAcrossElectrodes = cell(numRows,numCols,numProtocols);
psthNormalizedMeanAcrossElectrodes = cell(numRows,numCols,numProtocols);

if isfield(params,'extraEpochsFR')
    numExtraEpochs = length(params.extraEpochsFR);
    FRPerTrialEpoch = cell(numRows,numCols,numProtocols,numExtraEpochs);
    FRAcrossTrialsEpoch = cell(numRows,numCols,numProtocols,numExtraEpochs);
    semFRAcrossTrialsEpoch = cell(numRows,numCols,numProtocols,numExtraEpochs);
    meanFRAcrossTrialsElectrodesEpoch = cell(numRows,numCols,numProtocols,numExtraEpochs);
    semFRAcrossTrialsElectrodesEpoch = cell(numRows,numCols,numProtocols,numExtraEpochs);
end

% get the baseline and stimulus indices
% Fs = round(1/(params.timeVals(2)-params.timeVals(1)));
% BLRange = uint16((params.FRblmax-params.FRblmin)*Fs);
% blPos = find(params.timeVals>=params.FRblmin,1)+ (1:BLRange);
% STRange = uint16((params.FRstmax-params.FRstmin)*Fs);
% stPos = find(params.timeVals>=params.FRstmin,1)+ (1:STRange);

blPos = [params.FRblmin params.FRblmax];
stPos = [params.FRstmin params.FRstmax];

% Compute FR and related measures
for i = 1:numProtocols
    for r = 1:numRows
        for c = 1:numCols
            disp([' FR >> processing: protocol ' num2str(i) '/' num2str(numProtocols) ', row ' num2str(c) '/' num2str(numCols) ', col ' num2str(r) '/' num2str(numRows)]);
            if computepsth
                [psth{r,c,i},params.FRxs] = cellfun(@(x) getPSTH(x,params.FRtimebin,[params.timeVals(1) params.timeVals(end)],params.FRsmooothSigmaMS), data{r,c,i}, 'uniformoutput', 0);
            end
            FRPerTrialBL{r,c,i} = cellfun(@(x) getSpikeRate(x,blPos), data{r,c,i}, 'uniformoutput',0);
            FRPerTrialST{r,c,i} = cellfun(@(x) getSpikeRate(x,stPos), data{r,c,i}, 'uniformoutput',0);
            
            FRAcrossTrialsBL{r,c,i} = cellfun(@(x) nanmean(x), FRPerTrialBL{r,c,i}, 'uniformoutput',0);
            FRAcrossTrialsST{r,c,i} = cellfun(@(x) nanmean(x), FRPerTrialST{r,c,i}, 'uniformoutput',0);
            semFRAcrossTrialsBL{r,c,i} = cellfun(@(x) nanstd(x)/sqrt(length(x)), FRPerTrialBL{r,c,i}, 'uniformoutput',0);
            semFRAcrossTrialsST{r,c,i} = cellfun(@(x) nanstd(x)/sqrt(length(x)), FRPerTrialST{r,c,i}, 'uniformoutput',0);
            
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common')
                % weighted mean: the number of repeats per electrode may 
                % vary. In such a case one may take a weighted mean 
                % This approach is required only when you have unequal
                % repeats per electrode i.e. when the badTrials selection
                % strategy is not the 'common' badTrials strategy
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:)));
                meanFRAcrossTrialsElectrodesBL{r,c,i} = cell2mat(FRAcrossTrialsBL{r,c,i})*weights./sum(weights);
                meanFRAcrossTrialsElectrodesST{r,c,i} = cell2mat(FRAcrossTrialsST{r,c,i})*weights./sum(weights);
                
                % ##TODO Not sure how to take weighted std. Have kept it
                % the same as non-weighted std
                semFRAcrossTrialsElectrodesBL{r,c,i} = nanstd(cell2mat(FRAcrossTrialsBL{r,c,i}))/sqrt(numElecs);
                semFRAcrossTrialsElectrodesST{r,c,i} = nanstd(cell2mat(FRAcrossTrialsST{r,c,i}))/sqrt(numElecs);
                
                psthMeanAcrossElectrodes{r,c,i} = weights.*cell2mat(psth(r,c,i))./sum(weights);
                psthNormalizedMeanAcrossElectrodes{r,c,i} = weights.*cell2mat(psthNormalized(r,c,i))./sum(weights);
            else
                meanFRAcrossTrialsElectrodesBL{r,c,i} = nanmean(cell2mat(FRAcrossTrialsBL{r,c,i}));
                meanFRAcrossTrialsElectrodesST{r,c,i} = nanmean(cell2mat(FRAcrossTrialsST{r,c,i}));
                semFRAcrossTrialsElectrodesBL{r,c,i} = nanstd(cell2mat(FRAcrossTrialsBL{r,c,i}))/sqrt(numElecs);
                semFRAcrossTrialsElectrodesST{r,c,i} = nanstd(cell2mat(FRAcrossTrialsST{r,c,i}))/sqrt(numElecs);
                
                psthMeanAcrossElectrodes{r,c,i} = nanmean(cell2mat(psth(r,c,i)),1);
                psthNormalizedMeanAcrossElectrodes{r,c,i} = nanmean(cell2mat(psthNormalized(r,c,i)),1);
            end
            
            %**************************************************************
            % compute for extra epochs if they have been passed
            if isfield(params,'extraEpochsFR')
                for k=1:numExtraEpochs
                    timePos = [params.extraEpochsFR{k}(1) params.extraEpochsFR{k}(2)];
                    FRPerTrialEpoch{r,c,i,k} = cellfun(@(x) getSpikeRate(x,timePos), data{r,c,i}, 'uniformoutput',0);
                    FRAcrossTrialsEpoch{r,c,i,k} = cellfun(@(x) nanmean(x), FRPerTrialEpoch{r,c,i,k}, 'uniformoutput',0);
                    semFRAcrossTrialsEpoch{r,c,i,k} = cellfun(@(x) nanstd(x)/sqrt(length(x)), FRAcrossTrialsEpoch{r,c,i,k}, 'uniformoutput',0);

                    if params.weightedmean && ~strcmpi(params.badtrialsoption,'common')
                        % weighted mean: the number of repeats per electrode may 
                        % vary. In such a case one may take a weighted mean 
                        % This approach is required only when you have unequal
                        % repeats per electrode i.e. when the badTrials selection
                        % strategy is not the 'common' badTrials strategy
                        weights = cell2mat(squeeze(params.numRepeats(r,c,i,:)));
                        meanFRAcrossTrialsElectrodesEpoch{r,c,i,k} = cell2mat(FRAcrossTrialsEpoch{r,c,i,k})*weights./sum(weights);

                        % ##TODO Not sure how to take weighted std. Have kept it
                        % the same as non-weighted std
                        semFRAcrossTrialsElectrodesEpoch{r,c,i,k} = nanstd(cell2mat(FRAcrossTrialsEpoch{r,c,i,k}))/sqrt(numElecs);
                    else
                        meanFRAcrossTrialsElectrodesEpoch{r,c,i,k} = nanmean(cell2mat(FRAcrossTrialsEpoch{r,c,i,k}));
                        semFRAcrossTrialsElectrodesEpoch{r,c,i,k} = nanstd(cell2mat(FRAcrossTrialsEpoch{r,c,i,k}))/sqrt(numElecs);
                    end
                end
            end 
        end
    end
end

% normalized psth calculation
% get normalized PSTH : every electrode's psth is divided by
% its max psth value across all conditions
lengthpsth = size(psth{1},2);

for i = 1:numProtocols
    numElecs = size(data{1,1,i},2);
    psthmatrix = reshape([psth{:,:,i}],numElecs,lengthpsth,numRows,numCols);
    maxpsth = max(max(max(psthmatrix,[],2),[],3),[],4); % max across time, row and col i.e max for each electrode across all time and conditions
    for r = 1:numRows
        for c = 1:numCols
            psthNormalized(r,c,i) = cellfun(@(x) x.*repmat(1./maxpsth,1,lengthpsth), psth(r,c,i), 'uniformoutput',0); 
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common')
                % weighted mean: the number of repeats per electrode may 
                % vary. In such a case one may take a weighted mean 
                % This approach is required only when you have unequal
                % repeats per electrode i.e. when the badTrials selection
                % strategy is not the 'common' badTrials strategy
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:)));
                psthNormalizedMeanAcrossElectrodes{r,c,i} = weights.*cell2mat(psthNormalized(r,c,i))./sum(weights);
            else
                psthNormalizedMeanAcrossElectrodes{r,c,i} = nanmean(cell2mat(psthNormalized(r,c,i)),1);
            end
        end
    end
end

% Compute common baseline across conditions per electrode
commonBaselineFRAcrossTrials = cell(1,numProtocols);
for i = 1:numProtocols
    numElecs = size(data{1,1,i},2);
    commonBaselineFRAcrossTrials{i} = mean(reshape(cell2mat([FRAcrossTrialsBL{:,:,i}]),numElecs,numRows*numCols),2);
end

% Compute firing rate modulation
dFRAcrossTrials = cell(numRows, numCols, numProtocols);
dFRAcrossTrialsCommonBL = cell(numRows, numCols, numProtocols);
modFRAcrossTrials = cell(numRows, numCols, numProtocols);
modFRAcrossTrialsCommonBL = cell(numRows, numCols, numProtocols);
for i = 1:numProtocols
    for r = 1:numRows
        for c = 1:numCols
            dFRAcrossTrials{r,c,i} = cell2mat(FRAcrossTrialsST{r,c,i})-cell2mat(FRAcrossTrialsBL{r,c,i});
            dFRAcrossTrialsCommonBL{r,c,i} = cell2mat(FRAcrossTrialsST{r,c,i})-commonBaselineFRAcrossTrials{i}';
            
            modFRAcrossTrials{r,c,i} = cell2mat(cellfun(@(x,y) (x-y)./(x+y), FRAcrossTrialsST{r,c,i}, FRAcrossTrialsBL{r,c,i}, 'uniformoutput',0));
            modFRAcrossTrialsCommonBL{r,c,i} = cell2mat(cellfun(@(x,y) (x-y)./(x+y), FRAcrossTrialsST{r,c,i}, num2cell(commonBaselineFRAcrossTrials{i}'), 'uniformoutput',0));
        end
    end
end

%==========================================================================
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFRAcrossTrialsElectrodesSessionsBL = zeros(numRows,numCols);
meanFRAcrossTrialsElectrodesSessionsST = zeros(numRows,numCols);
dFRAcrossTrialsSessions = zeros(numRows,numCols);
dFRAcrossTrialsSessionsCommonBL = zeros(numRows,numCols);
modFRAcrossTrialsSessions = zeros(numRows,numCols);
modFRAcrossTrialsSessionsCommonBL = zeros(numRows,numCols);
for r = 1:numRows
    for c = 1:numCols
        if numProtocols>1 % take average across protocols when analyzing multiple protocols
            if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                clear weights
                weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:))),2)';
                numElecsSession = cell2mat(squeeze(cellfun(@(x) size(x,2), dFRAcrossTrials(r,c,:),'uniformoutput',0)));

                meanFRAcrossTrialsElectrodesSessionsBL(r,c,:) = squeeze(weights*(squeeze(cell2mat(meanFRAcrossTrialsElectrodesBL(r,c,:))))/sum(weights));
                meanFRAcrossTrialsElectrodesSessionsST(r,c,:) = squeeze(weights*(squeeze(cell2mat(meanFRAcrossTrialsElectrodesST(r,c,:))))/sum(weights));

                dFRAcrossTrialsSessions(r,c,:) = sum(cell2mat((cellfun(@(x,y) x.*y./(weights*numElecsSession),squeeze(dFRAcrossTrials(r,c,:)),num2cell(weights)','uniformoutput',0))'));
                dFRAcrossTrialsSessionsCommonBL(r,c,:) = sum(cell2mat((cellfun(@(x,y) x.*y./(weights*numElecsSession),squeeze(dFRAcrossTrialsCommonBL(r,c,:)),num2cell(weights)','uniformoutput',0))'));
                
                modFRAcrossTrialsSessions(r,c,:) = sum(cell2mat((cellfun(@(x,y) x.*y./(weights*numElecsSession),squeeze(modFRAcrossTrials(r,c,:)),num2cell(weights)','uniformoutput',0))'));
                modFRAcrossTrialsSessionsCommonBL(r,c,:) = sum(cell2mat((cellfun(@(x,y) x.*y./(weights*numElecsSession),squeeze(modFRAcrossTrialsCommonBL(r,c,:)),num2cell(weights)','uniformoutput',0))'));
            else
                meanFRAcrossTrialsElectrodesSessionsBL(r,c,:) = squeeze(nanmean(squeeze(cell2mat(meanFRAcrossTrialsElectrodesBL(r,c,:))),1)); % squeeze(cell2mat(.)) is a single column matrix, hence mean across 1st dimension 
                meanFRAcrossTrialsElectrodesSessionsST(r,c,:) = squeeze(nanmean(squeeze(cell2mat(meanFRAcrossTrialsElectrodesST(r,c,:))),1)); % squeeze(cell2mat(.)) is a single column matrix, hence mean across 1st dimension 

                dFRAcrossTrialsSessions(r,c,:) = squeeze(nanmean([dFRAcrossTrials{r,c,:}],2));
                dFRAcrossTrialsSessionsCommonBL(r,c,:) = squeeze(nanmean([dFRAcrossTrialsCommonBL{r,c,:}],2));
                
                modFRAcrossTrialsSessions(r,c,:) = squeeze(nanmean([modFRAcrossTrials{r,c,:}],2));
                modFRAcrossTrialsSessionsCommonBL(r,c,:) = squeeze(nanmean([modFRAcrossTrialsCommonBL{r,c,:}],2));

            end
        else
            meanFRAcrossTrialsElectrodesSessionsBL(r,c,:) = meanFRAcrossTrialsElectrodesBL{r,c,:};
            meanFRAcrossTrialsElectrodesSessionsST(r,c,:) = meanFRAcrossTrialsElectrodesST{r,c,:};

            dFRAcrossTrialsSessions(r,c,:) = nanmean(dFRAcrossTrials{r,c,:},2);
            dFRAcrossTrialsSessionsCommonBL(r,c,:) = nanmean(dFRAcrossTrialsCommonBL{r,c,:},2);
            
            modFRAcrossTrialsSessions(r,c,:) = nanmean(modFRAcrossTrials{r,c,:},2);
            modFRAcrossTrialsSessionsCommonBL(r,c,:) = nanmean(modFRAcrossTrialsCommonBL{r,c,:},2);
        end
        
        if isfield(params,'extraEpochsFR')
            meanFRAcrossTrialsElectrodesSessionsEpoch = zeros(numRows,numCols,numExtraEpochs);
            for k=1:numExtraEpochs
                if numProtocols>1
                    if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                        clear weights
                        weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:))),2)';
                        meanFRAcrossTrialsElectrodesSessionsEpoch(r,c,k,:) = squeeze(weights*(squeeze(cell2mat(meanFRAcrossTrialsElectrodesEpoch(r,c,:,k))))/sum(weights));
                    else
                        meanFRAcrossTrialsElectrodesSessionsEpoch(r,c,k,:) = squeeze(nanmean(squeeze(cell2mat(meanFRAcrossTrialsElectrodesEpoch(r,c,:,k))),1)); % squeeze(cell2mat(.)) is a single column matrix, hence mean across 1st dimension 
                    end
                else
                    meanFRAcrossTrialsElectrodesSessionsEpoch(r,c,k,:) = meanFRAcrossTrialsElectrodesEpoch{r,c,:,k};
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
% Plot FRs
% Raw FRs: plot baseline and stimulus period FRs
% Check if plot handles have been passed, else create a new figure and plot
% handles
if ~isfield(params,'hFRFig') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hFRFig = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

if isfield(params,'hFRFig')
    psthMeanAcrossElectrodesSessions = zeros(numRows,numCols,length(params.FRxs));
    psthNormalizedMeanAcrossElectrodesSessions = zeros(numRows,numCols,length(params.FRxs));
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:))),2)';

                    psthMeanAcrossElectrodesSessions(r,c,:) = squeeze((squeeze(cell2mat(psthMeanAcrossElectrodes(r,c,:)))*weights')./sum(weights));
                    psthNormalizedMeanAcrossElectrodesSessions(r,c,:) = squeeze((squeeze(cell2mat(psthNormalizedMeanAcrossElectrodes(r,c,:)))*weights')./sum(weights));

                else
                    psthMeanAcrossElectrodesSessions(r,c,:) = squeeze(nanmean(squeeze(cell2mat(psthMeanAcrossElectrodes(r,c,:))),2));
                    psthNormalizedMeanAcrossElectrodesSessions(r,c,:) = squeeze(nanmean(squeeze(cell2mat(psthNormalizedMeanAcrossElectrodes(r,c,:))),2));

                end
            else
                psthMeanAcrossElectrodesSessions(r,c,:) = psthMeanAcrossElectrodes{r,c,:};
                psthNormalizedMeanAcrossElectrodesSessions(r,c,:) = psthNormalizedMeanAcrossElectrodes{r,c,:};
            end
            
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hFRFig(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hFRFig(x,y));
            end
            
            if strcmpi(params.psthType,'raw')
                plot(params.FRxs,squeeze(psthMeanAcrossElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
            else % normalized
                plot(params.FRxs,squeeze(psthNormalizedMeanAcrossElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
            end
            set(gca,'nextplot','add');
%             setAxesProperties(params.hFRFig(x,y),params,x,y,rows);
%             showTitleAndN(params.hFRFig(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end
    set(params.hFRFig,'ylim',ylimits);
    set(params.hFRFig,'xlim',[params.FRxs(1) params.FRxs(end)]);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean FR across trials, electrodes and sessions','FitBoxToText','on');
    end    

    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hFRFig(r,c));
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
        if ~isfield(params,'hFRMarginalsFigRow') && newfiguresflag
            figure;
            plotsPos = [0.08 0.15 0.87 0.8];
            params.hFRMarginalsFigRow = getPlotHandles(numRows,1,plotsPos);
        end
        
        if isfield(params,'hFRMarginalsFigRow')
            for r = 1:numRows % one plot per row
                subplot(params.hFRMarginalsFigRow(r,1));
                for c = 1:numCols
                    if strcmpi(params.psthType,'raw')
                        plot(params.FRxs,squeeze(psthMeanAcrossElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
                    else % normalized
                        plot(params.FRxs,squeeze(psthNormalizedMeanAcrossElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
                    end
                    set(gca,'nextplot','add');
                    title(params.titleStringY{r,1});
                end

                if r==numRows
                    setAxesProperties(gca,params,1,1,1);
                end
                if r==ceil(numRows/2)
                    legend(params.titleStringX);
                    title(params.titleStringY{r,1});
                end
            end
            if newfiguresflag
                annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
                annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean FR across trials, electrodes and sessions','FitBoxToText','on');
            end
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % column-wise
        if ~isfield(params,'hFRMarginalsFigCol') && newfiguresflag
            figure;
            plotsPos = [0.08 0.15 0.87 0.8];
            params.hFRMarginalsFigCol = getPlotHandles(1,numCols,plotsPos);
        end
        if isfield(params,'hFRMarginalsFigCol')
            for c = 1:numCols % one plot per col
                subplot(params.hFRMarginalsFigCol(1,c));
                for r = 1:numRows
                    if strcmpi(params.psthType,'raw')
                        plot(params.FRxs,squeeze(psthMeanAcrossElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
                    else % normalized
                        plot(params.FRxs,squeeze(psthNormalizedMeanAcrossElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
                    end
                    set(gca,'nextplot','add');
                    title(params.titleStringX{1,c});
                end
                if c==1
                    setAxesProperties(gca,params,1,1,1);
                end
                if c==ceil(numCols/2)
                    legend(params.titleStringY{:,c});
                end
            end
            if newfiguresflag
                annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
                annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean FR across trials, electrodes and sessions','FitBoxToText','on');
            end
        end
    end
end
%==========================================================================
%--------------------------------------------------------------------------
% Plot Trends for the FR measures (band power, peak frequency)
%--------------------------------------------------------------------------

if isfield(params,'showtrends') && isfield(params,'polar')
    
    if params.showtrends
        if ~isfield(params,'hdFRTrend') && newfiguresflag
            figure;
            plotsPos = [0.08 0.15 0.87 0.8];
            params.hdFRTrend = getPlotHandles(numRows,1,plotsPos);
        end
        
        if isfield(params,'hdFRTrend')
            % Plot change in power in the selected band
            subplot(params.hdFRTrend); 
            for r = 1:numRows
                for c = 1:numCols
                    if strcmpi(params.FRbaselinetype,'conditionwise')
                        plot(params.valsUnique{params.paramx}(c),mean(squeeze(dFRAcrossTrialsSessions(r,c,:))),'o','color',params.colorsList((r-1)*numCols+c,:));
                    else
                        plot(params.valsUnique{params.paramx}(c),mean(squeeze(dFRAcrossTrialsSessionsCommonBL(r,c,:))),'o','color',params.colorsList((r-1)*numCols+c,:));
                    end
                    set(gca,'nextplot','add');
                end
                if strcmpi(params.FRbaselinetype,'conditionwise')
                    plot(params.valsUnique{params.paramx},mean(squeeze(dFRAcrossTrialsSessions(r,:,:)),2),'color',[0.8 0.8 0.8]/r);
                else
                    plot(params.valsUnique{params.paramx},mean(squeeze(dFRAcrossTrialsSessionsCommonBL(r,:,:)),2),'color',[0.8 0.8 0.8]/r);
                end
                linename = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
                text(0.1,0.95-0.1*r,linename,'fontsize',6,'color',[0.8 0.8 0.8]/r,'unit','normalized');
            end
            title('Change in firing rate (sp/s)');
            ylabel('dFR');
            line([0 params.valsUnique{params.paramx}(numCols)],[0 0],'color',[0.2 0.7 0.4]);
            setAxesBasicProperties(gca);
        end
        
        if ~isfield(params,'hmodFRTrend') && newfiguresflag
            figure;
            plotsPos = [0.08 0.15 0.87 0.8];
            params.hmodFRTrend = getPlotHandles(numRows,1,plotsPos);
        end
        
        if isfield(params,'hmodFRTrend')
            % Plot peak frequency in the selected band
            subplot(params.hmodFRTrend);
            for r = 1:numRows
                for c = 1:numCols
                    if strcmpi(params.FRbaselinetype,'conditionwise')
                        plot(params.valsUnique{params.paramx}(c),nanmean(squeeze(modFRAcrossTrialsSessions(r,c,:))),'o','color',params.colorsList((r-1)*numCols+c,:));
                    else
                        plot(params.valsUnique{params.paramx}(c),nanmean(squeeze(modFRAcrossTrialsSessionsCommonBL(r,c,:))),'o','color',params.colorsList((r-1)*numCols+c,:));
                    end
                    set(gca,'nextplot','add');
                end
                if strcmpi(params.FRbaselinetype,'conditionwise')
                    plot(params.valsUnique{params.paramx},nanmean(squeeze(modFRAcrossTrialsSessions(r,:,:)),2),'color',[0.8 0.8 0.8]/r);
                else
                    plot(params.valsUnique{params.paramx},nanmean(squeeze(modFRAcrossTrialsSessionsCommonBL(r,:,:)),2),'color',[0.8 0.8 0.8]/r);
                end
                linename = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
                text(0.1,0.95-0.1*r,linename,'fontsize',6,'color',[0.8 0.8 0.8]/r,'unit','normalized');
            end
            title('Modulation of firing rate (sp/s)');
            ylabel('modFR');
            line([0 params.valsUnique{params.paramx}(numCols)],[0 0],'color',[0.2 0.7 0.4]);
            setAxesBasicProperties(gca);
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
                if ~isfield(params,'hdFRPolar') && newfiguresflag
                    handlesDontExist = 1;
                else
                    handlesDontExist = 0;
                end
                if handlesDontExist
                    figure;
                    params.hdFRPolar = subplot(2,2,1);
                end
                
                if strcmpi(params.FRbaselinetype,'conditionwise')
                    filldata = mean(squeeze(dFRAcrossTrialsSessions(r,:,:)),2);
                else
                    filldata = mean(squeeze(dFRAcrossTrialsSessionsCommonBL(r,:,:)),2);
                end
                numTheta = numOris-extraDataPoints;
                hValsList = params.valsUnique{5}(1:numTheta); % ori = hue
                if extraDataPoints
                    extraPts = filldata(end-extraDataPoints:end);
                    makePolarPlotHSV(params.hdFRPolar,filldata',numTheta,hValsList,'linear','extra',num2cell(extraPts),'numCircles',numCircles,'overwrite',0,'colorsList',params.colorsList,'palettelevel',numCircles+1);
                else
                    numCircles = ceil(max(max(filldata)))-floor(min(min(filldata)));
                    makePolarPlotHSV(params.hdFRPolar,filldata',numTheta,hValsList,'linear','numCircles',numCircles,'overwrite',0,'colorsList',params.colorsList,'palettelevel',numCircles+1);
                end
                title('change in FR');
                setAxesBasicProperties(params.hdFRPolar);
                set(params.hdFRPolar,'yticklabel',[]);
                set(params.hdFRPolar,'xticklabel',[]);

                if handlesDontExist
                    polarfigname = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
                    set(gcf,'numbertitle','off','name',polarfigname);
                end
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if ~isfield(params,'hmoddFRPolar') && newfiguresflag
                    handlesDontExist = 1;
                else
                    handlesDontExist = 0;
                end
                if handlesDontExist
                    figure;
                    params.hmoddFRPolar = subplot(2,2,2);
                end
                
                if strcmpi(params.FRbaselinetype,'conditionwise')
                    filldata = nanmean(squeeze(modFRAcrossTrialsSessions(r,:,:)),2);
                else
                    filldata = nanmean(squeeze(modFRAcrossTrialsSessionsCommonBL(r,:,:)),2);
                end
                numTheta = numOris-extraDataPoints;
                hValsList = params.valsUnique{5}(1:numTheta); % ori = hue
                if extraDataPoints
                    extraPts = filldata(end-extraDataPoints:end);
                    makePolarPlotHSV(params.hmoddFRPolar,filldata',numTheta,hValsList,'linear','extra',num2cell(extraPts),'numCircles',numCircles,'overwrite',1,'colorsList',params.colorsList,'palettelevel',numCircles+1);
                else
                    numCircles = ceil(max(max(filldata)))-floor(min(min(filldata)));
                    makePolarPlotHSV(params.hmoddFRPolar,filldata',numTheta,hValsList,'linear','numCircles',numCircles,'overwrite',1,'colorsList',params.colorsList,'palettelevel',numCircles+1);
                end
                title('modulation of FR');
                setAxesBasicProperties(params.hmoddFRPolar);
                set(params.hmoddFRPolar,'yticklabel',[]);
                set(params.hmoddFRPolar,'xticklabel',[]);

                if handlesDontExist
                    polarfigname = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
                    set(gcf,'numbertitle','off','name',polarfigname);
                end
            end

        else % rows are represented as hues
            for c=1:numCases
                if ~isfield(params,'hdFRPolar') && newfiguresflag
                    handlesDontExist = 1;
                else
                    handlesDontExist = 0;
                end
                if handlesDontExist
                    figure;
                    params.hdFRPolar = subplot(2,2,1);
                end
                
                if strcmpi(params.FRbaselinetype,'conditionwise')
                    filldata = mean(squeeze(dFRAcrossTrialsSessions(:,c,:)),2);
                else
                    filldata = mean(squeeze(dFRAcrossTrialsSessionsCommonBL(:,c,:)),2);
                end
                numTheta = numOris-extraDataPoints;
                hValsList = params.valsUnique{5}(1:numTheta); % ori = hue
                if extraDataPoints
                    extraPts = filldata(end-extraDataPoints:end);
                    makePolarPlotHSV(params.hdFRPolar,filldata',numTheta,hValsList,'linear','extra',num2cell(extraPts),'numCircles',numCircles,'overwrite',0,'colorsList',params.colorsList,'palettelevel',numCircles+1);
                else
                    numCircles = ceil(max(max(filldata)))-floor(min(min(filldata)));
                    makePolarPlotHSV(params.hdFRPolar,filldata',numTheta,hValsList,'linear','numCircles',numCircles,'overwrite',0,'colorsList',params.colorsList,'palettelevel',numCircles+1);
                end
                title('change in FR');
                setAxesBasicProperties(params.hdFRPolar);
                set(params.hdFRPolar,'yticklabel',[]);
                set(params.hdFRPolar,'xticklabel',[]);

                if handlesDontExist
                    polarfigname = [params.paramsString{params.paramx} ':' num2str(params.valsUnique{params.paramx}(c))];
                    set(gcf,'numbertitle','off','name',polarfigname);
                end
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if ~isfield(params,'hmoddFRPolar') && newfiguresflag
                    handlesDontExist = 1;
                else
                    handlesDontExist = 0;
                end
                if handlesDontExist
                    figure;
                    params.hmoddFRPolar = subplot(2,2,2);
                end
                
                if strcmpi(params.FRbaselinetype,'conditionwise')
                    filldata = nanmean(squeeze(modFRAcrossTrialsSessions(:,c,:)),2);
                else
                    filldata = nanmean(squeeze(modFRAcrossTrialsSessionsCommonBL(:,c,:)),2);
                end
                numTheta = numOris-extraDataPoints;
                hValsList = params.valsUnique{5}(1:numTheta); % ori = hue
                if extraDataPoints
                    extraPts = filldata(end-extraDataPoints:end);
                    makePolarPlotHSV(params.hmoddFRPolar,filldata',numTheta,hValsList,'linear','extra',num2cell(extraPts),'numCircles',numCircles,'overwrite',1,'colorsList',params.colorsList,'palettelevel',numCircles+1);
                else
                    numCircles = ceil(max(max(filldata)))-floor(min(min(filldata)));
                    makePolarPlotHSV(params.hmoddFRPolar,filldata',numTheta,hValsList,'linear','numCircles',numCircles,'overwrite',1,'colorsList',params.colorsList,'palettelevel',numCircles+1);
                end
                title('modulation of FR');
                setAxesBasicProperties(params.hmoddFRPolar);
                set(params.hmoddFRPolar,'yticklabel',[]);
                set(params.hmoddFRPolar,'xticklabel',[]);

                if handlesDontExist
                    polarfigname = [params.paramsString{params.paramx} ':' num2str(params.valsUnique{params.paramx}(c))];
                    set(gcf,'numbertitle','off','name',polarfigname);
                end
            end
        end
    end
end

%==========================================================================
%--------------------------------------------------------------------------
% Make topoplots
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
end
%==========================================================================
%==========================================================================
% return the specified fr data
%--------------------------------------------------------------------------
if returnData
    frdata.psth = psth;
    frdata.psthNormalized = psthNormalized;
    frdata.psthMeanAcrossElectrodes = psthMeanAcrossElectrodes;
    frdata.psthNormalizedMeanAcrossElectrodes = psthNormalizedMeanAcrossElectrodes;
    frdata.frPerTrialBL = FRPerTrialBL;
    frdata.frPerTrialST = FRPerTrialST;
    frdata.frBL = FRAcrossTrialsBL;
    frdata.frST = FRAcrossTrialsST;
    frdata.frMeanBL = meanFRAcrossTrialsElectrodesBL;
    frdata.frMeanST = meanFRAcrossTrialsElectrodesST;
    if strcmpi(params.baselinestrategy,'conditionwise')
        frdata.diffFR = dFRAcrossTrials;
        frdata.frModulation = modFRAcrossTrials;
    elseif strcmpi(params.baselinestrategy,'common')        
        frdata.diffFR = dFRAcrossTrialsCommonBL;
        frdata.frModulation = modFRAcrossTrialsCommonBL;
        frdata.commonBaselineFRAcrossTrials = commonBaselineFRAcrossTrials;
    end
    
    if isfield(params,'extraEpochsFR')
        frdata.frPerTrialEpoch = FRPerTrialEpoch;
        frdata.frEpoch = FRAcrossTrialsEpoch;
        frdata.frMeanEpoch = meanFRAcrossTrialsElectrodesEpoch;
    end
    
    if drawplots && newfiguresflag
        frdata.psthMeanAcrossElectrodesSessions = psthMeanAcrossElectrodesSessions;
        frdata.psthNormalizedMeanAcrossElectrodesSessions = psthNormalizedMeanAcrossElectrodesSessions;
        frdata.frMeanSessionsBL = meanFRAcrossTrialsElectrodesSessionsBL;
        frdata.frMeanSessionsST = meanFRAcrossTrialsElectrodesSessionsST;
        
        if strcmpi(params.baselinestrategy,'conditionwise')
            frdata.diffFRMeanSessions = dFRAcrossTrialsSessions;
            frdata.frModulationMeanSessions = modFRAcrossTrialsSessions;
        elseif strcmpi(params.baselinestrategy,'common')        
            frdata.diffFRMeanSessions = dFRAcrossTrialsSessionsCommonBL;
            frdata.frModulationMeanSessions = modFRAcrossTrialsSessionsCommonBL;
        end

        if isfield(params,'extraEpochsFR')
            frdata.frMeanSessionsEpoch = meanFRAcrossTrialsElectrodesSessionsEpoch;
        end
    end
    frdata.FRxs = params.FRxs;
    frdata.params = params;
else
    frdata = [];
end


end