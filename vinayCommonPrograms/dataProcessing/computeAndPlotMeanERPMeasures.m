% plot mean ERP across electrodes, sessions, conditions
% Vinay Shirhatti, 28 September 2016
% 
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
% *************************************************************************

function [params,erpdata] = computeAndPlotMeanERPMeasures(data,params,drawplots,newfiguresflag,returnData)

if ~exist('drawplots','var'); drawplots=1; end
if ~exist('newfiguresflag','var'); newfiguresflag=0; end
if ~exist('returnData','var'); returnData=1; end

%--------------------------------------------------------------------------
% Read all sizes
numRows = size(data,1);
numCols = size(data,2);
numProtocols = size(data,3);
numElecs = size(data{1},2);

baselineCorrection = 1; % if you want DC corrected ERP

if ~isfield(params,'maxPlotsAlongX')
    params.maxPlotsAlongX = 6; % this variable stores the maximum number of plots to be drawn along the horizontal axis
end

% assign timeVals if not present
if ~isfield(params,'timeVals')
    params.timeVals = 1:2048;
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
% Compute and plot ERPs
%--------------------------------------------------------------------------
% Initialize variables
ERPAcrossTrials = cell(numRows,numCols,numProtocols);
stdERPAcrossTrials = cell(numRows,numCols,numProtocols);
meanERPAcrossTrialsElectrodes = cell(numRows,numCols,numProtocols);
stimPeriodMin = cell(numRows,numCols,numProtocols);
stimPeriodMax = cell(numRows,numCols,numProtocols);
rmsChangeBLtoST = cell(numRows,numCols,numProtocols);

% get indices for the erp epoch used for calculating measures
Fs = round(1/(params.timeVals(2)-params.timeVals(1)));
BLRange = uint16((params.blmax-params.blmin)*Fs);
blPos = find(params.timeVals>=params.blmin,1)+ (1:BLRange);
ERPRange = uint16((params.erpmax-params.erpmin)*Fs);
erpPos = find(params.timeVals>=params.erpmin,1)+ (1:ERPRange);

% Compute ERP and related measures
for i = 1:numProtocols
    for r = 1:numRows
        for c = 1:numCols
            disp([' processing: protocol ' num2str(i) '/' num2str(numProtocols) ', row ' num2str(c) '/' num2str(numCols) ', col ' num2str(r) '/' num2str(numRows)]);
            ERPAcrossTrials{r,c,i} = cell2mat(cellfun(@mean,data{r,c,i}','UniformOutput',0)); % ERP for every electrode averaged across the trials
            if baselineCorrection % subtract the mean baseline (i.e baseline DC) from ERP
                blDC = mean(ERPAcrossTrials{r,c,i}(:,blPos),2);
                ERPAcrossTrials{r,c,i} = ERPAcrossTrials{r,c,i} - repmat(blDC,1,size(ERPAcrossTrials{r,c,i},2));
            end
            stdERPAcrossTrials{r,c,i} = cell2mat(cellfun(@(x) std(x,[],1),data{r,c,i}','UniformOutput',0)); % do we need baseline correction for std cal? I think we will get the same value even then
            
            % Mean across electrodes
            if params.weightedmean && ~strcmpi(params.badtrialsoption,'common') 
                % weighted mean: the number of repeats per electrode may 
                % vary. In such a case one may take a weighted mean 
                % This approach is required only when you have unequal
                % repeats per electrode i.e. when the badTrials selection
                % strategy is not the 'common' badTrials strategy
                weights = cell2mat(squeeze(params.numRepeats(r,c,i,:)))';
                meanERPAcrossTrialsElectrodes{r,c,i} = weights*ERPAcrossTrials{r,c,i}/sum(weights); % this is a matrix multiplication
                % weights dim: 1 x numElecs => contains the number of 
                % repeats for each electrode for this condition and session
                % ERPAcrossTrials{r,c,i} dim: numElecs x numTimePoints =>
                % contains the mean across trials for each electrode.
                % Through this matrix operation we weight each erp by the
                % number of repeats used to compute it and then take the
                % mean of all erps across electrodes. The resultant matrix
                % is of dimensions 1 x numTimePoints i.e. the weighted
                % average of erps across electrodes
            else
                meanERPAcrossTrialsElectrodes{r,c,i} = mean(ERPAcrossTrials{r,c,i},1);
            end
            
            %_________________ERP Measures___________________
            stimPeriodMin{r,c,i} = min(ERPAcrossTrials{r,c,i}(:,erpPos)'); % min in the erp period
            stimPeriodMax{r,c,i} = max(ERPAcrossTrials{r,c,i}(:,erpPos)'); % max in the erp period
            rmsChangeBLtoST{r,c,i} = rms(ERPAcrossTrials{r,c,i}(:,erpPos)') - rms(ERPAcrossTrials{r,c,i}(:,blPos)'); % rms change from baseline to erp period
        end
    end
end


%==========================================================================
if drawplots
%==========================================================================
%--------------------------------------------------------------------------
% Plot ERPs
%--------------------------------------------------------------------------

% Check if plot handles have been passed, else create a new figure and plot
% handles
if ~isfield(params,'hERPFig') && newfiguresflag
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX % the number of plots allowed along the x axis are less than the total cases then adjust the number of rows and columns accordingly
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0; % no rearrangement otherwise
    end
    params.hERPFig = getPlotHandles(rows,cols,plotsPos);
end

if isfield(params,'hERPFig')
    numTimePoints = length(params.timeVals);
    meanERPAcrossTrialsElectrodesSessions = zeros(numRows,numCols,numTimePoints);
    
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols>1 % take average across protocols when analyzing multiple protocols
                if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = sum(cell2mat(squeeze(params.numRepeats(r,c,:,:)))')';
                    meanERPAcrossTrialsElectrodesSessions(r,c,:) = squeeze((squeeze(cell2mat(meanERPAcrossTrialsElectrodes(r,c,:)))*weights)/sum(weights));
                else
                    meanERPAcrossTrialsElectrodesSessions(r,c,:) = squeeze(nanmean(squeeze(cell2mat(meanERPAcrossTrialsElectrodes(r,c,:))),2));
                end
            else
                meanERPAcrossTrialsElectrodesSessions(r,c,:) = meanERPAcrossTrialsElectrodes{r,c,:};
            end
            if ~subplotsRearrangement
                x = r; y = c;
                subplot(params.hERPFig(x,y));
            else
                x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                y = mod(c,cols);
                if y==0; y = cols; end
                subplot(params.hERPFig(x,y));
            end
            plot(params.timeVals,squeeze(meanERPAcrossTrialsElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
            setAxesProperties(params.hERPFig(x,y),params,x,y,rows);
            showTitleAndN(params.hERPFig(x,y),params,r,c);

            if y==1 && x==rows
                ylimits = ylim;
            end
            drawnow;
        end
    end
    set(params.hERPFig,'ylim',ylimits);
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean ERP across trials, electrodes and sessions','FitBoxToText','on');
    end
    
    % delete unused axes
    for r=1:rows
        for c = 1:cols
            if cols*(r-1)+c>numRows*numCols
                axisHandle = findall(params.hERPFig(r,c));
                delete(axisHandle);
            end
        end
    end 
end

%==========================================================================
%--------------------------------------------------------------------------
% Plot marginals
%--------------------------------------------------------------------------
if params.showmarginals
    % Row-wise
    if ~isfield(params,'hMarginalsFigRow') && newfiguresflag
        figure;
        plotsPos = [0.08 0.15 0.87 0.8];
        params.hMarginalsFigRow = getPlotHandles(numRows,1,plotsPos);
    end
    
    if isfield(params,'hMarginalsFigRow')
        for r = 1:numRows % one plot per row
            subplot(params.hMarginalsFigRow(r,1));
            for c = 1:numCols
                plot(params.timeVals,squeeze(meanERPAcrossTrialsElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
                set(gca,'nextplot','add');
                title(params.titleStringY{r,1});
            end

            if r==numRows
                setAxesProperties(gca,params,1,1,1);
            end
            if r==ceil(numRows/2)
                legend(gca,params.titleStringY{:,c});
                title(params.titleStringY{r,1});
            end
        end
    end
    
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean ERP across trials, electrodes and sessions','FitBoxToText','on');
    end
    
    % column-wise
    if ~isfield(params,'hMarginalsFigCol') && newfiguresflag
        figure;
        plotsPos = [0.08 0.15 0.87 0.8];
        params.hMarginalsFigCol = getPlotHandles(1,numCols,plotsPos);
    end
    
    if isfield(params,'hMarginalsFigCol')
        for c = 1:numCols % one plot per col
            subplot(params.hMarginalsFigCol(1,c));
            for r = 1:numRows
                plot(params.timeVals,squeeze(meanERPAcrossTrialsElectrodesSessions(r,c,:)),'color',params.colorsList((r-1)*numCols+c,:));
                set(gca,'nextplot','add');
                title(params.titleStringX{1,c});
            end
            if c==1
                setAxesProperties(gca,params,1,1,1);
            end
            if c==ceil(numCols/2)
                legend(gca,params.titleStringY{:,c});
            end
        end
    end
    
    if newfiguresflag
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean ERP across trials, electrodes and sessions','FitBoxToText','on');
    end
end

%==========================================================================
%--------------------------------------------------------------------------
% Plot Trends for the ERP measures (min, max, rms)
%--------------------------------------------------------------------------

meanstimPeriodMin = zeros(numRows,numCols);
meanstimPeriodMax = zeros(numRows,numCols);
meanrmsChangeBLtoST = zeros(numRows,numCols);
if params.showtrends || params.polar % these computations are done and plots are plotted when any one of 'trends' or 'polar' is selected (or both are selected)
    for r = 1:numRows
        for c = 1:numCols
            if numProtocols==1 % take average across protocols only when analyzing multiple protocols
                meanstimPeriodMin(r,c,:) = stimPeriodMin{r,c,:};
                meanstimPeriodMax(r,c,:) = stimPeriodMax{r,c,:};
                meanrmsChangeBLtoST(r,c,:) = rmsChangeBLtoST{r,c,:};
            else
                 if params.weightedmean % take weighted average across sessions. The weights are the total number of repeats for that given session
                    clear weights
                    weights = cell2mat(squeeze(params.numRepeats(r,c,:,:)))'; % contains the repeats per electrode per session
                    % In this case the value for each electrode is to be
                    % weighted by the repeats for that session and then the 
                    % mean value has to be calculated for each electrode across
                    % the sessions. For eg.:
                    % stimPeriodMin is numElec x numProtocols
                    % weights is also numElec x numProtocols as calculated here
                    % so we take elementwise product so that each value is
                    % weighted by the number of repeats used to calculate that
                    % value. Then we sum across the sessions and divide by the
                    % total number of repeats across the sessions

                    meanstimPeriodMin(r,c,:) = nanmean(sum((squeeze(cell2mat(stimPeriodMin(r,c,:))).*weights),2)./sum(weights')');
                    meanstimPeriodMax(r,c,:) = nanmean(sum((squeeze(cell2mat(stimPeriodMax(r,c,:))).*weights),2)./sum(weights')');
                    meanrmsChangeBLtoST(r,c,:) = nanmean(sum((squeeze(cell2mat(rmsChangeBLtoST(r,c,:))).*weights),2)./sum(weights')');
                 else
                    meanstimPeriodMin(r,c,:) = squeeze(mean(nanmean(squeeze(cell2mat(stimPeriodMin(r,c,:))),2)));
                    meanstimPeriodMax(r,c,:) = squeeze(mean(nanmean(squeeze(cell2mat(stimPeriodMax(r,c,:))),2)));
                    meanrmsChangeBLtoST(r,c,:) = squeeze(mean(nanmean(squeeze(cell2mat(rmsChangeBLtoST(r,c,:))),2)));
                 end
            end
        end
    end
    
    if newfiguresflag
        figure;
        % Plot ERP min
        subplot(3,1,1); 
        for r = 1:numRows
            for c = 1:numCols
                plot(params.valsUnique{params.paramx}(c),squeeze(meanstimPeriodMin(r,c,:)),'o','color',params.colorsList((r-1)*numCols+c,:));
                set(gca,'nextplot','add');
            end
            plot(params.valsUnique{params.paramx},squeeze(meanstimPeriodMin(r,:,:)),'color',[0.8 0.8 0.8]/r);
            linename = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
            text(0.1,0.95-0.1*r,linename,'fontsize',6,'color',[0.8 0.8 0.8]/r,'unit','normalized');
        end
        title(['Minimum amplitude during stimulus period [' num2str(params.erpmin) ' ' num2str(params.erpmax) '] s']);
        ylabel('uV');
        setAxesBasicProperties(gca);

        % Plot ERP max
        subplot(3,1,2);
        for r = 1:numRows
            for c = 1:numCols
                plot(params.valsUnique{params.paramx}(c),squeeze(meanstimPeriodMax(r,c,:)),'o','color',params.colorsList((r-1)*numCols+c,:));
                set(gca,'nextplot','add');
            end
            plot(params.valsUnique{params.paramx},squeeze(meanstimPeriodMax(r,:,:)),'color',[0.8 0.8 0.8]/r);
            linename = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
            text(0.1,0.95-0.1*r,linename,'fontsize',6,'color',[0.8 0.8 0.8]/r,'unit','normalized');
        end
        title(['Maximum amplitude during stimulus period [' num2str(params.erpmin) ' ' num2str(params.erpmax) '] s']);
        ylabel('uV');
        setAxesBasicProperties(gca);

        % Plot ERP rms change
        subplot(3,1,3);
        for r = 1:numRows
            for c = 1:numCols
                plot(params.valsUnique{params.paramx}(c),squeeze(meanrmsChangeBLtoST(r,c,:)),'o','color',params.colorsList((r-1)*numCols+c,:));
                set(gca,'nextplot','add');
            end
            plot(params.valsUnique{params.paramx},squeeze(meanrmsChangeBLtoST(r,:,:)),'color',[0.8 0.8 0.8]/r);
            linename = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
            text(0.1,0.95-0.1*r,linename,'fontsize',6,'color',[0.8 0.8 0.8]/r,'unit','normalized');
        end
        title(['RMS change from baseline [' num2str(params.blmin) ' ' num2str(params.blmax) '] s to stimulus period [' num2str(params.erpmin) ' ' num2str(params.erpmax) '] s']);
        ylabel('uV^2');
        setAxesBasicProperties(gca);

        annotation('textbox',[0.4 0.97 0.4 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
    end
end

%--------------------------------------------------------------------------
% Make polar plots
%--------------------------------------------------------------------------
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
            if ~isfield(params,'hERPMinPolar')
                handlesDontExist = 1;
            else
                handlesDontExist = 0;
            end
            if handlesDontExist
                figure;
                params.hERPMinPolar = subplot(2,3,1);
                params.hERPMaxPolar = subplot(2,3,2);
                params.hERPRMSPolar = subplot(2,3,3);
            end
            polarfigname = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
            filldata = squeeze(meanstimPeriodMin(r,:,:));
            numTheta = numOris-extraDataPoints;
            hValsList = params.valsUnique{5}(1:numTheta); % ori = hue
            makePolarPlotHSV(params.hERPMinPolar,filldata,numTheta,hValsList,'abs','extra',extraDataPoints);
            title('ERP Min');
            setAxesBasicProperties(params.hERPMinPolar);
            set(params.hERPMinPolar,'yticklabel',[]);
            set(params.hERPMinPolar,'xticklabel',[]);

            filldata = squeeze(meanstimPeriodMax(r,:,:));
            makePolarPlotHSV(params.hERPMaxPolar,filldata,numTheta,hValsList,'abs','extra',extraDataPoints);
            title('ERP Max');
            setAxesBasicProperties(params.hERPMaxPolar);
            set(params.hERPMaxPolar,'yticklabel',[]);
            set(params.hERPMaxPolar,'xticklabel',[]);

            filldata = squeeze(meanrmsChangeBLtoST(r,:,:));
            makePolarPlotHSV(params.hERPRMSPolar,filldata,numTheta,hValsList,'abs','extra',extraDataPoints);
            title('ERP RMS');
            setAxesBasicProperties(params.hERPRMSPolar);
            set(params.hERPRMSPolar,'yticklabel',[]);
            set(params.hERPRMSPolar,'xticklabel',[]);

            set(gcf,'numbertitle','off','name',polarfigname);
        end
        
    else % rows are represented as hues
        for c=1:numCases
            if ~isfield(params,'hERPMinPolar')
                handlesDontExist = 1;
            else
                handlesDontExist = 0;
            end
            if handlesDontExist
                figure;
                params.hERPMinPolar = subplot(2,3,1);
                params.hERPMaxPolar = subplot(2,3,2);
                params.hERPRMSPolar = subplot(2,3,3);
            end
            
            polarfigname = [params.paramsString{params.paramx} ':' num2str(params.valsUnique{params.paramx}(c))];
            filldata = (squeeze(meanstimPeriodMin(:,c,:)))';
            numTheta = numOris-extraDataPoints;
            hValsList = params.valsUnique{5}(1:numTheta); % ori = hue
            makePolarPlotHSV(params.hERPMinPolar,filldata,numTheta,hValsList,'abs','extra',extraDataPoints);
            title('ERP Min');
            setAxesBasicProperties(params.hERPMinPolar);
            set(params.hERPMinPolar,'yticklabel',[]);
            set(params.hERPMinPolar,'xticklabel',[]);

            filldata = (squeeze(meanstimPeriodMax(:,c,:)))';
            makePolarPlotHSV(params.hERPMaxPolar,filldata,numTheta,hValsList,'abs','extra',extraDataPoints);
            title('ERP Max');
            setAxesBasicProperties(params.hERPMaxPolar);
            set(params.hERPMaxPolar,'yticklabel',[]);
            set(params.hERPMaxPolar,'xticklabel',[]);

            filldata = (squeeze(meanrmsChangeBLtoST(:,c,:)))';
            makePolarPlotHSV(params.hERPRMSPolar,filldata,numTheta,hValsList,'abs','extra',extraDataPoints);
            title('ERP RMS');
            setAxesBasicProperties(params.hERPRMSPolar);
            set(params.hERPRMSPolar,'yticklabel',[]);
            set(params.hERPRMSPolar,'xticklabel',[]);

            set(gcf,'numbertitle','off','name',polarfigname);
        end
    end

end

%==========================================================================
%--------------------------------------------------------------------------
% Make topoplots
%--------------------------------------------------------------------------
if params.topo
    for r=1:numRows
        if numProtocols>1 % take average across protocols only when analyzing multiple protocols
            if params.weightedmean && strcmpi(params.badtrialsoption,'common')
                clear weights
                weights = cell2mat(squeeze(params.numRepeats(r,:,:,:))); % contains the repeats per electrode per session
                % weights: numCols x numProtocols
                for c=1:numCols
                    stimPeriodMinPerElectrodeAcrossSessions(r,:,c) = weights(c,:)*cell2mat(squeeze(stimPeriodMin(r,c,:)))./sum(weights(c,:));
                    stimPeriodMaxPerElectrodeAcrossSessions(r,:,c) = weights(c,:)*cell2mat(squeeze(stimPeriodMax(r,c,:)))./sum(weights(c,:));
                    rmsChangeBLtoSTPerElectrodeAcrossSessions(r,:,c) = weights(c,:)*cell2mat(squeeze(rmsChangeBLtoST(r,c,:)))./sum(weights(c,:));
                end
            elseif params.weightedmean && ~strcmpi(params.badtrialsoption,'common') % individual bad trials
                clear weights
                for c=1:numCols
                    weights = cell2mat(squeeze(params.numRepeats(r,c,:,:)));
                    % weight: numProtocols x numElecs 
                    % We use sum here because it is not a matrix multiplication
                    stimPeriodMinPerElectrodeAcrossSessions(r,:,c) = sum(weights.*cell2mat(squeeze(stimPeriodMin(r,c,:))),1)./sum(weights,1);
                    stimPeriodMaxPerElectrodeAcrossSessions(r,:,c) = sum(weights.*cell2mat(squeeze(stimPeriodMax(r,c,:))),1)./sum(weights,1);
                    rmsChangeBLtoSTPerElectrodeAcrossSessions(r,:,c) = sum(weights.*cell2mat(squeeze(rmsChangeBLtoST(r,c,:))),1)./sum(weights,1);
                end
            else
                stimPeriodMinPerElectrodeAcrossSessions(r,:,:) = reshape((nanmean(squeeze(cell2mat(stimPeriodMin(r,:,:))),2)),numElecs,numCols);
                stimPeriodMaxPerElectrodeAcrossSessions(r,:,:) = reshape((nanmean(squeeze(cell2mat(stimPeriodMax(r,:,:))),2)),numElecs,numCols);
                rmsChangeBLtoSTPerElectrodeAcrossSessions(r,:,:) = reshape((nanmean(squeeze(cell2mat(rmsChangeBLtoST(r,:,:))),2)),numElecs,numCols);
            end
            
        else
            for c=1:numCols
                stimPeriodMinPerElectrodeAcrossSessions(r,:,c) = stimPeriodMin{r,c,:};
                stimPeriodMaxPerElectrodeAcrossSessions(r,:,c) = stimPeriodMax{r,c,:};
                rmsChangeBLtoSTPerElectrodeAcrossSessions(r,:,c) = rmsChangeBLtoST{r,c,:};
            end
        end
    end
    
    % decide the number of topoplots and assign data accordingly
    clear dataAllConditions
    if strcmpi(params.erptype,'all')
        dataAllConditions{1} = rmsChangeBLtoSTPerElectrodeAcrossSessions;
        dataAllConditions{2} = stimPeriodMinPerElectrodeAcrossSessions;
        dataAllConditions{3} = stimPeriodMaxPerElectrodeAcrossSessions;
        figname{1} = 'topoplot_rmsERP';
        figname{2} = 'topoplot_minERP';
        figname{3} = 'topoplot_maxERP';
    elseif strcmpi(params.erptype,'rms')
        dataAllConditions{1} = rmsChangeBLtoSTPerElectrodeAcrossSessions;
        figname{1} = 'topoplot_rmsERP';
    elseif strcmpi(params.erptype,'min')
        dataAllConditions{1} = stimPeriodMinPerElectrodeAcrossSessions;
        figname{1} = 'topoplot_minERP';
    elseif strcmpi(params.erptype,'max')
        dataAllConditions{1} = stimPeriodMaxPerElectrodeAcrossSessions;
        figname{1} = 'topoplot_maxERP';
    end
    
    % Draw the topoplots
    handlesDontExist = 0;
    for k = 1:length(dataAllConditions)  
        if ~isfield(params,'hTopoFig')
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
            params.hTopoFig{k} = getPlotHandles(rows,cols,plotsPos);
        end
        
        if ~iscell(params.hTopoFig)
            tempArray = params.hTopoFig;
            params = rmfield(params,'hTopoFig');
            params.hTopoFig{1} = tempArray;
        end
        
        colormap(jet);

        for r = 1:numRows
            for c = 1:numCols
                if ~subplotsRearrangement
                    x = r; y = c;
                    subplot(params.hTopoFig{k}(x,y));
                else
                    x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                    y = mod(c,cols);
                    if y==0; y = cols; end
                    subplot(params.hTopoFig{k}(x,y));
                end

                datatoplot = zeros(1,length(params.chanlocs));
                datatoplot(:,params.electrodes) = squeeze(dataAllConditions{k}(r,:,c));
                topoplot(datatoplot,params.chanlocs,'electrodes','off','drawaxis','off',...
                    'plotchans',params.electrodes,'hcolor',[0.25 0.25 0.25]);
                whitebg([0 0 0]);
                showTitleAndN(params.hTopoFig{k}(x,y),params,r,c,0,-0.25,0.94); % show title and n outside the top left part of topoplot
                drawnow;
            end
        end
        set(params.hTopoFig{k},'clim',[min(min(min(dataAllConditions{k}))) max(max(max(dataAllConditions{k})))]);
        annotation('textbox',[0.1 0.97 0.6 0.02],'String',[params.reftype{1} ' reference, Electrodes: ' num2str(params.electrodes)],'FitBoxToText','on');
        annotation('textbox',[0.75 0.97 0.23 0.02],'String','Mean across trials and sessions','FitBoxToText','on');
        colorbar('position',[plotsPos(1)+plotsPos(3)+0.005 plotsPos(2) 0.015 plotsPos(4)]);
        set(gcf,'numbertitle','off','name',figname{k});
        
        % delete unused axes
        for r=1:rows
            for c = 1:cols
                if cols*(r-1)+c>numRows*numCols
                    axisHandle = findall(params.hTopoFig{k}(r,c));
                    delete(axisHandle);
                end
            end
        end 
    end    
end

%==========================================================================
%--------------------------------------------------------------------------
% Show the good stimuli at every electrode, for every condition
%--------------------------------------------------------------------------

if params.showgoodstim
    for en = 1:numElecs
        figure;
        plotsPos = [0.08 0.15 0.87 0.8];
        if numCols>params.maxPlotsAlongX
            cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
            subplotsRearrangement = 1;
        else
            cols = numCols; rows = numRows; subplotsRearrangement = 0;
        end
        params.hElecGoodStimFig{en} = getPlotHandles(rows,cols,plotsPos);
        for r = 1:numRows
            for c = 1:numCols
                if ~subplotsRearrangement
                    x = r; y = c;
                    subplot(params.hElecGoodStimFig{en}(x,y));
                else
                    x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
                    y = mod(c,cols);
                    if y==0; y = cols; end
                    subplot(params.hElecGoodStimFig{en}(x,y));
                end
                
                plot(params.timeVals,data{r,c,i}{en});
                axis('tight');
            end
        end
        
        % delete unused axes
        for r=1:rows
            for c = 1:cols
                if cols*(r-1)+c>numRows*numCols
                    axisHandle = findall(params.hElecGoodStimFig{en}(r,c));
                    delete(axisHandle);
                end
            end
        end
        set(gcf,'numbertitle','off','name',['electrode ' num2str(params.electrodes(en))]);
    end
end
%--------------------------------------------------------------------------
end
%==========================================================================

if returnData
    erpdata.ERPAcrossTrials = ERPAcrossTrials;
    erpdata.stdERPAcrossTrials = stdERPAcrossTrials;
    erpdata.meanERPAcrossTrialsElectrodes = meanERPAcrossTrialsElectrodes;
    erpdata.stimPeriodMin = stimPeriodMin;
    erpdata.stimPeriodMax = stimPeriodMax;
    erpdata.rmsChangeBLtoST = rmsChangeBLtoST;
    erpdata.timeVals = params.timeVals;
    erpdata.params = params;
else
    erpdata = [];
end
%==========================================================================

end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%==========================================================================
% Additional Functions
%==========================================================================

function setAxesProperties(hAxes,params,row,col,rows)

if ~isfield(params,'setAxesProperties')
    params.setAxesProperties = 1;
end


if ~isfield(params,'cleanAxes')
    params.cleanAxes = 0;
end

if params.setAxesProperties
    set(hAxes,'xlim',[params.tmin params.tmax]);
    set(hAxes,'tickLength',[0.03 0.01]);
    set(hAxes,'fontsize',20,'fontweight','bold');
    set(hAxes,'xtick',[0 0.8],'tickdir','out');
    set(hAxes,'tickdir','out');
    set(hAxes,'box','off');

    if col==1 && row==rows
        ylabel(hAxes,'ERP (uV)');
        xlabel(hAxes,'Time (s)');
    else
        set(hAxes,'yticklabel',[]);
        set(hAxes,'xticklabel',[]);
    end
elseif params.cleanAxes
    set(hAxes,'yticklabel',[]);
    set(hAxes,'xticklabel',[]);
else
    setAxesBasicProperties(hAxes);
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

function setAxesBasicProperties(hAxes)
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