% make plots for EEG data
% Vinay Shirhatti, 26 September 2016
%
% modified for plotting LFP data
% 22 October 2016
% *************************************************************************

function [params,dataOut] = makeplotsLFP(params)

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% determine the electrodes to analyze (LFP)
useAsAnalogInput = 0; % initialize flag to use channels as ainp instead of elecs
for i=1:length(params.electrodeSelectionLFP)
    switch params.electrodeSelectionLFP{i}
        case {'none'} % no electrode selected => no analysis required on electrode signals. 
            % This choice can be taken while doing only stimulus related operations, eg. drawing the different stimuli in the protocol 
            params.electrodes{i} = [];
            params.tf = 0; params.psd = 0; params.erp = 0;
        case {'all'} % select all electrodes for analyses
            params.electrodes{i} = 1:96;
        case {'highRMS'}
            rfDataFile = [params.subjectName{1} params.gridType 'RFData.mat'];
            load(rfDataFile);
            params.electrodes{i} = highRMSElectrodes;
        case {'lowRMS'}
            rfDataFile = [params.subjectName{1} params.gridType 'RFData.mat'];
            load(rfDataFile);
            params.electrodes{i} = lowRMSElectrodes;
        case {'highLowRMS'}
            rfDataFile = [params.subjectName{1} params.gridType 'RFData.mat'];
            load(rfDataFile);
            params.electrodes{i} = [highRMSElectrodes lowRMSElectrodes];
        case {'manual'} % chosen electrodes are stored in params.electrodes
            params.electrodes{i} = params.electrodesManualLFP{i};
        case {'ainp'}
            if isfield(params,'electrodesManualLFP')
                if ~isempty(params.electrodesManualLFP{i})
                    params.electrodes{i} = params.electrodesManualLFP{i};
                else
                    params.electrodes{i} = 1:6;
                end
            else
                params.electrodes{i} = 1:6;
            end
            useAsAnalogInput = 1;
        otherwise
            params.electrodes{i} = 1:96;
    end
    params.electrodes{i} = setdiff(params.electrodes{i},params.ignoreElectrodes{i});
    params.allelectrodes{i} = params.electrodes{i};
end

% determine the electrodes to analyze (FR)
for i=1:length(params.electrodeSelectionSpikes)
    switch params.electrodeSelectionSpikes{i}
        case {'none'} % no electrode selected => no analysis required on electrode signals. 
            % This choice can be taken while doing only stimulus related operations, eg. drawing the different stimuli in the protocol 
            params.electrodesFR{i} = [];
            params.tf = 0; params.psd = 0; params.erp = 0;
        case {'all'} % select all electrodes for analyses
            params.electrodesFR{i} = 1:96;
        case {'highRMS'}
            rfDataFile = [params.subjectName{1} params.gridType 'RFData.mat'];
            load(rfDataFile);
            params.electrodesFR{i} = highRMSElectrodes;
        case {'lowRMS'}
            rfDataFile = [params.subjectName{1} params.gridType 'RFData.mat'];
            load(rfDataFile);
            params.electrodesFR{i} = lowRMSElectrodes;
        case {'highLowRMS'}
            rfDataFile = [params.subjectName{1} params.gridType 'RFData.mat'];
            load(rfDataFile);
            params.electrodesFR{i} = [highRMSElectrodes lowRMSElectrodes];
        case {'manual'} % chosen electrodes are stored in params.electrodes
            params.electrodesFR{i} = params.electrodesManualSpikes{i};
        otherwise
            params.electrodesFR{i} = 1:96;
    end
    params.electrodesFR{i} = setdiff(params.electrodesFR{i},params.ignoreElectrodes{i});
    params.allelectrodesFR{i} = params.electrodesFR{i};
end

if ~isfield(params,'rejectBadElectrodes')
    params.rejectBadElectrodes = 1; % by default reject bad electrodes
end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% parameters to be varied/pooled along the axes and parameters with fixed
% values.
% Create a cell matrix of size rows x cols where rows = number of values of
% paramy and cols = number of values of x. For every position we create a
% cell matrix of 1x9 to store the values of 9 parameters to be considered
% for that position. The goodPos for that condition are then determined 
% based on these parameter values.
% The set of parameter values to be considered for analyses are stored in 2
% variables: params.valsUnique (which has all the single values) and 
% params.valsUniquePool (which has the pooling sets). The total number of
% cases is the combination of single value cases and the pooling sets.
x = params.valsUnique{params.paramx}; % individual values along x
y = params.valsUnique{params.paramy}; % individual values along y
if sum(cell2mat(params.boolPool)) % at least one parameter has pooling sets
    xp = params.valsUniquePool{params.paramx}; % pooling sets along x
    yp = params.valsUniquePool{params.paramy}; % pooling sets along y
else
    xp = []; yp = [];
end
numRows = length(y)+size(yp,1);
numCols = length(x)+size(xp,1);
valsMatrix = cell(numRows,numCols,9); % initialize
for c = 1:numCols
    for r = 1:numRows
        for i=1:9 % assign values for 9 parameters for each condition
            if i==params.paramx % parameter along x i.e. columns
                if c <= length(x)
                    valsMatrix{r,c,i} = params.valsUnique{i}(c);
                else % add the pooling sets
                    valsMatrix{r,c,i} = params.valsUniquePool{i}(c-length(x),:);
                end
            elseif i==params.paramy % parameter along y i.e. rows
                if r <= length(y)
                    valsMatrix{r,c,i} = params.valsUnique{i}(r);
                else % add the pooling sets
                    valsMatrix{r,c,i} = params.valsUniquePool{i}(r-length(y),:);
                end
            else % other parameters
                valsMatrix{r,c,i} = params.valsUnique{i}; % other parameters are fixed across all conditions
            end
            valsMatrix{r,c,i}(isnan(valsMatrix{r,c,i})) = []; % remove NaNs
        end
    end
end

% make the legend strings
params.titleStringX = cell(numRows,numCols);
params.titleStringY = cell(numRows,numCols);
for c = 1:numCols
    for r = 1:numRows
        for i=1:9 % check 9 parameters for each condition
            params.titleStringX{r,c} = [params.paramsString{params.paramx} ':' num2str(params.valsUnique{params.paramx}(c))];
            params.titleStringY{r,c} = [params.paramsString{params.paramy} ':' num2str(params.valsUnique{params.paramy}(r))];
        end
    end
end

% Read the badTrials strategy and set flags accordingly
if strcmpi(params.badtrialsoption{1},'common')
    removeCommonBadtrials=1; % to remove the common badTrials only
    removeIndividualBadtrials = 0;
elseif strcmpi(params.badtrialsoption{1},'individual')
    removeCommonBadtrials=0; 
    removeIndividualBadtrials = 1; % to remove individual badTrials
else % none
    removeCommonBadtrials=0; 
    removeIndividualBadtrials = 0;
end

numProtocols = length(params.subjectslist);
allBadTrials = cell(1,numProtocols); % stores the list of badTrials for every electrode, every protocol
params.badTrials = cell(1,numProtocols); % stores the list of common badTrials for every protocol
elecData = cell(numRows,numCols,numProtocols); % the MAIN DATA cell array. Stores the data for every electrode for every condition and protocol
spikeData = cell(numRows,numCols,numProtocols); % the MAIN SPIKE DATA array. Stores the data for every electrode for every condition and protocol
psthVals = cell(numRows,numCols,numProtocols); % PSTH values from the spike data
params.spikeChannelNumber = cell(1,numProtocols); % spike channel numbers (this list could be different from the params.electrodes)
params.unitIDList = cell(1,numProtocols); % list of spike unit IDs

if ~isfield(params,'saveDataFolder')
    saveDataFolder = mkdir(fullfile(params.folderSourceString,'data','savedData'));
else
    saveDataFolder = params.saveDataFolder;
end
tag = [params.subjectslist(:) params.expdateslist(:) params.protocolslist(:)]';
tag = [tag{:}];
if ~isfield(params,'tagUser')
    tagUser = [];
else
    tagUser = params.tagUser;
end
tag = [tag tagUser];
tagLFP = cell(1,numProtocols);
tagFR = cell(1,numProtocols);

for i=1:numProtocols
    if strcmpi(params.electrodeSelectionLFP{i},'manual')
        if length(params.electrodes{i})==1
            tagLFP{i} = [tag '_goodPosElecData_' params.electrodeSelectionLFP{i} '_elec' num2str(params.electrodes{i}(1))];
        else
            tagLFP{i} = [tag '_goodPosElecData_' params.electrodeSelectionLFP{i} '_' num2str(length(params.electrodes{i})) 'elecs'];
        end
    elseif strcmpi(params.electrodeSelectionLFP{i},'ainp')
        if length(params.electrodes{i})==1
            tagLFP{i} = [tag '_goodPosElecData_' params.electrodeSelectionLFP{i} num2str(params.electrodes{i}(1))];
        else
            tagLFP{i} = [tag '_goodPosElecData_' params.electrodeSelectionLFP{i} '_' num2str(length(params.electrodes{i})) 'ainps'];
        end
    else
        tagLFP{i} = [tag '_goodPosElecData_' params.electrodeSelectionLFP{i} '_' num2str(length(params.electrodes{i})) 'elecs'];
    end
    if strcmpi(params.electrodeSelectionSpikes{i},'manual')
        if length(params.electrodesFR{i})==1
            tagFR{i} = [tag '_goodPosFR_' params.electrodeSelectionSpikes{i} '_elec' num2str(params.electrodesFR{i}(1))];
        else
            tagFR{i} = [tag '_goodPosFR_' params.electrodeSelectionSpikes{i} '_' num2str(length(params.electrodesFR{i})) 'elecs'];
        end
    else
        tagFR{i} = [tag '_goodPosFR_' params.electrodeSelectionSpikes{i} '_' num2str(length(params.electrodesFR{i})) 'elecs'];
    end
end
tagLFP = [tagLFP{:}];
tagFR = [tagFR{:}];

processElecData = 1;
processFRData = 1;
if ~(params.erp || params.psd || params.tf)
    processElecData = 0;
end
if ~params.fr; processFRData = 0; end

savedDataExists = 0;
savedDataExistsFR = 0;
if params.loadSavedData && ~params.overwriteElecData && processElecData       
    savedfile = fullfile(saveDataFolder,[tagLFP '.mat']);
    if exist(savedfile,'file')
        disp('Loading saved electrode data >>>>>');
        savedElecData = load(savedfile,'elecData','params');
        elecData = savedElecData.elecData;
        params.timeVals = savedElecData.params.timeVals;
        params.analyzeElectrodes = savedElecData.params.analyzeElectrodes;
        params.badTrials = savedElecData.params.badTrials;
        params.badElecs = savedElecData.params.badElecs;
        params.badImpedance = savedElecData.params.badImpedance;
        params.numRepeats = savedElecData.params.numRepeats;
        if removeIndividualBadtrials
            params.elecbadtrials = savedElecData.params.elecbadtrials;
        end
        savedDataExists = 1;
        clear savedElecData
    end
end
if params.loadSavedData && ~params.overwriteFRData && processFRData
    savedfile = fullfile(saveDataFolder,[tagFR '.mat']);
    if exist(savedfile,'file')
        disp('Loading saved FR data >>>>>');
        savedFRData = load(savedfile,'spikeData','psthVals','params');
        spikeData = savedFRData.spikeData;
        psthVals = savedFRData.psthVals;
        params.FRxs = savedFRData.params.FRxs;
        params.spikeChannelNumber = savedFRData.params.spikeChannelNumber;
        params.unitIDList = savedFRData.params.unitIDList;
        params.analyzeElectrodesFR = savedFRData.params.analyzeElectrodesFR;
        params.badTrials = savedFRData.params.badTrials;
        params.badElecs = savedFRData.params.badElecs;
        params.badImpedance = savedFRData.params.badImpedance;
        params.numRepeatsFR = savedFRData.params.numRepeatsFR;
        params.FRtimebin = savedFRData.params.FRtimebin;
        params.snrThreshold = savedFRData.params.snrThreshold;
        params.psthThreshold = savedFRData.params.psthThreshold;
        params.removeUnitID255 = savedFRData.params.removeUnitID255;
        params.useAsMultiunit = savedFRData.params.useAsMultiunit;
        if removeIndividualBadtrials
            params.elecbadtrialsFR = savedFRData.params.elecbadtrialsFR;
        end
        savedDataExistsFR = 1;
        clear savedFRData
    end 
end

for i = 1:numProtocols
    if ~savedDataExists && processElecData
        if ~isempty(params.electrodes{i}) && (params.erp || params.psd || params.tf)
            disp(['Processing protocol number: ' num2str(i) ' out of ' num2str(length(params.subjectslist))]);

            numElecs = length(params.allelectrodes{i});
            if strcmpi(params.badtrialsoption{1},'common')
                goodPos = cell(numRows,numCols,numProtocols); % initialize goodPos
            elseif strcmpi(params.badtrialsoption{1},'individual')
                goodPos = cell(numRows,numCols,numProtocols,numElecs); % initialize goodPos
            else % none
                goodPos = cell(numRows,numCols,numProtocols); % initialize goodPos
            end

            subject = params.subjectslist{i}; % read subject
            expdate = params.expdateslist{i}; % read expdate
            protocol = params.protocolslist{i}; % read protocol

            folderName = fullfile(params.folderSourceString,'data',subject,params.gridType,expdate,protocol);
            folderSegment = fullfile(folderName,'segmentedData');
            load(fullfile(folderSegment,'LFP','lfpInfo.mat')); % load lfpInfo for timeVals

            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Remove all bad electrodes
            % params.badElecs : electrodes with statistically high number of
            %                   badTrials
            % params.badImpedance : electrodes with impedance outside the
            %                       acceptable range
            % params.badCrosstalkElectrodes : electrodes with a high ('red')
            %                                 level of crosstalk
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            [allBadTrials{i}, params.badTrials{i}, params.badElecs{i}] = loadBadTrials(fullfile(folderSegment,'badTrials.mat')); % load the list of bad trials for that protocol

            if strcmpi(params.electrodeSelectionLFP{i},'ainp')
                params.analyzeElectrodes{i} = params.electrodes{i};
            else
                if ~params.rejectBadElectrodes
                    disp('_______NOT REMOVING ANY ELECTRODE_________');
                    params.analyzeElectrodes{i} = params.electrodes{i};
                    params.badImpedance{i} = [];
                else
                    % remove bad impedance electrodes if any        
                    disp('**** removing bad impedance electrodes ****');
                    badElecs = getBadImpedanceElectrodes(params.folderSourceString,subject,expdate,params.gridType,params.impedanceCutoffs{i});
                    if ~isempty(badElecs)
                        disp(num2str(badElecs));
                    else
                        disp('No bad impedance electrodes');
                    end
                    params.badImpedance{i} = badElecs;

                    params.analyzeElectrodes{i} = setdiff(params.allelectrodes{i},params.badImpedance{i});
                    params.analyzeElectrodes{i} = setdiff(params.analyzeElectrodes{i},params.badElecs{i});

                    if isfield('params','badCrosstalkElectrodes')
                        params.analyzeElectrodes{i} = setdiff(params.analyzeElectrodes{i},params.badCrosstalkElectrodes{i});
                    end
                    params.electrodes{i} = params.analyzeElectrodes{i};
                end
            end
            disp(['***** ' num2str(length(params.analyzeElectrodes{i})) ' electrodes classified as good ******']);
            disp(num2str(params.analyzeElectrodes{i}));

%             if params.fr % Determine spiking units
%                 if ~isfield(params,'FRtimebin')
%                     params.FRtimebin = 10; % default: 10ms
%                 end
%                 if ~isfield(params,'snrThreshold')
%                     params.snrThreshold{i} = 0; % default: 0
%                 end
% 
%                 if ~isfield(params,'useSortedSpikesKMeans')
%                     [~,~,~,params.spikeChannelNumber{i},params.unitIDList{i}] = getMeanFiringRate(folderName,'electrodeList',params.electrodes{i},'goodPos',1,'timebin',params.FRtimebin,'snrThreshold',params.snrThreshold{i},'psthScreening','psthThreshold',params.psthThreshold{i},'useAllunitIDs','snrScreening');
%                 else
%                     if params.useSortedSpikesKMeans{i}
%                         [~,~,~,params.spikeChannelNumber{i},params.unitIDList{i}] = getMeanFiringRate(folderName,'electrodeList',params.electrodes{i},'goodPos',1,'timebin',params.FRtimebin,'snrThreshold',params.snrThreshold{i},'psthScreening','psthThreshold',params.psthThreshold{i},'useAllunitIDs','useSortedSpikesKMeans','snrScreening');
%                     else
%                         [~,~,~,params.spikeChannelNumber{i},params.unitIDList{i}] = getMeanFiringRate(folderName,'electrodeList',params.electrodes{i},'goodPos',1,'timebin',params.FRtimebin,'snrThreshold',params.snrThreshold{i},'psthScreening','psthThreshold',params.psthThreshold{i},'useAllunitIDs','snrScreening');
%                     end
%                 end
%             end

            % Now generate the data for each condition
            for c = 1:numCols
                for r = 1:numRows
                    disp(['***************** processing: protocol ' num2str(i) '/' num2str(numProtocols) ', col ' num2str(c) '/' num2str(numCols) ', row ' num2str(r) '/' num2str(numRows)]);

                    % get the goodPos for this particular condition
                    if strncmpi(protocol,'GRF',3)
                        goodPos{r,c,i} = getGoodPosGRF(folderName,squeeze(valsMatrix(r,c,:)),removeCommonBadtrials);
                    elseif strncmpi(protocol,'CRS',3)
                        goodPos{r,c,i} = getGoodPosCRS(folderName,squeeze(valsMatrix(r,c,:)),removeCommonBadtrials);
                    end
                    if isfield(params,'ignoreStimPos')
                        goodPos{r,c,i} = setdiff(goodPos{r,c,i},params.ignoreStimPos{i});
                    end

                    % if badTrials are considered individually for each electrode
                    % then compute the goodPos for each electrode
                    numElecs = length(params.electrodes{i});
                    if removeIndividualBadtrials
                        for en = 1:numElecs
                            if strcmpi(params.reftype{1},'single')
                                params.elecbadtrials{i,en} = allBadTrials{i}{params.electrodes{i}(en)};
                            else % bipolar - not worked out for LFP yet
                                params.elecbadtrials{i,en} = union(allBadTrials{i}{params.chanElecs(en,1)},allBadTrials{i}{params.chanElecs(en,2)});
                            end
                            goodPos{r,c,i,en} = setdiff(goodPos{r,c,i},params.elecbadtrials{i,en});
                            params.numRepeats{r,c,i,en} = size(goodPos{r,c,i,en},2); % store the number of repeats for this condition
                        end
                    else
                        params.numRepeats{r,c,i} = size(goodPos{r,c,i},2); % store the number of repeats for this condition
                    end

                    if params.erp || params.psd || params.tf
                        elecData{r,c,i} = getAnalogDataForElecAndPos(folderName,squeeze(goodPos(r,c,i,:)),params.electrodes{i},params.singleprecision,params.reftype,[],useAsAnalogInput);
                    end

%                     if params.fr
%                         if ~isfield(params,'useSortedSpikesKMeans')
%                             [spikeData{r,c,i},psthVals{r,c,i},params.FRxs,params.spikeChannelNumber{r,c,i},params.unitIDList{r,c,i}] = getMeanFiringRate(folderName,'electrodeList',params.spikeChannelNumber{i},'goodPos',squeeze(goodPos(r,c,i,:)),'timebin',params.FRtimebin,'unitIDList',params.unitIDList{i});
%                         else
%                             if params.useSortedSpikesKMeans{i}
%                                 [spikeData{r,c,i},psthVals{r,c,i},params.FRxs,params.spikeChannelNumber{r,c,i},params.unitIDList{r,c,i}] = getMeanFiringRate(folderName,'electrodeList',params.spikeChannelNumber{i},'goodPos',squeeze(goodPos(r,c,i,:)),'timebin',params.FRtimebin,'unitIDList',params.unitIDList{i},'useSortedSpikesKMeans');
%                             else
%                                 [spikeData{r,c,i},psthVals{r,c,i},params.FRxs,params.spikeChannelNumber{r,c,i},params.unitIDList{r,c,i}] = getMeanFiringRate(folderName,'electrodeList',params.spikeChannelNumber{i},'goodPos',squeeze(goodPos(r,c,i,:)),'timebin',params.FRtimebin,'unitIDList',params.unitIDList{i});
%                             end
%                         end
%                     end
                end
            end
            params.timeVals = timeVals;
        end
    end
end

if ~savedDataExists && params.saveDataFlag && processElecData % just load data if available and not generate the electrode data everytime
    disp('<><><><>~~ Saving goodPos analog data ~~<><><><> ');
    save(fullfile(saveDataFolder,tagLFP),'elecData','params')
%     numRepeats = params.numRepeats;
%     save('elecData.mat','elecData','numRepeats');
end

for i = 1:numProtocols
    if ~savedDataExistsFR && processFRData
        if ~isempty(params.electrodesFR{i}) && params.fr
            disp(['Processing FR, protocol number: ' num2str(i) ' out of ' num2str(length(params.subjectslist))]);

            numElecs = length(params.allelectrodesFR{i});
            if strcmpi(params.badtrialsoption{1},'common')
                goodPos = cell(numRows,numCols,numProtocols); % initialize goodPos
            elseif strcmpi(params.badtrialsoption{1},'individual')
                goodPos = cell(numRows,numCols,numProtocols,numElecs); % initialize goodPos
            else % none
                goodPos = cell(numRows,numCols,numProtocols); % initialize goodPos
            end

            subject = params.subjectslist{i}; % read subject
            expdate = params.expdateslist{i}; % read expdate
            protocol = params.protocolslist{i}; % read protocol

            folderName = fullfile(params.folderSourceString,'data',subject,params.gridType,expdate,protocol);
            folderSegment = fullfile(folderName,'segmentedData');
            load(fullfile(folderSegment,'LFP','lfpInfo.mat')); % load lfpInfo for timeVals

            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Remove all bad electrodes
            % params.badElecs : electrodes with statistically high number of
            %                   badTrials
            % params.badImpedance : electrodes with impedance outside the
            %                       acceptable range
            % params.badCrosstalkElectrodes : electrodes with a high ('red')
            %                                 level of crosstalk
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            [allBadTrials{i}, params.badTrials{i}, params.badElecs{i}] = loadBadTrials(fullfile(folderSegment,'badTrials.mat')); % load the list of bad trials for that protocol

            if ~params.rejectBadElectrodes
                disp('_______NOT REMOVING ANY ELECTRODE_________');
                params.analyzeElectrodes{i} = params.electrodes{i};
                params.badImpedance{i} = [];
            else
                % remove bad impedance electrodes if any        
                disp('**** removing bad impedance electrodes ****');
                badElecs = getBadImpedanceElectrodes(params.folderSourceString,subject,expdate,params.gridType,params.impedanceCutoffs{i});
                if ~isempty(badElecs)
                    disp(num2str(badElecs));
                else
                    disp('No bad impedance electrodes');
                end
                params.badImpedance{i} = badElecs;

                params.analyzeElectrodesFR{i} = setdiff(params.allelectrodesFR{i},params.badImpedance{i});
                params.analyzeElectrodesFR{i} = setdiff(params.analyzeElectrodesFR{i},params.badElecs{i});

                if isfield('params','badCrosstalkElectrodes')
                    params.analyzeElectrodesFR{i} = setdiff(params.analyzeElectrodesFR{i},params.badCrosstalkElectrodes{i});
                end
                params.electrodesFR{i} = params.analyzeElectrodesFR{i};
            end
            
            disp(['***** ' num2str(length(params.analyzeElectrodesFR{i})) ' electrodes classified as good ******']);
            disp(num2str(params.analyzeElectrodesFR{i}));

            if params.fr % Determine spiking units
                if ~isfield(params,'FRtimebin')
                    params.FRtimebin = 10; % default: 10ms
                end
                if ~isfield(params,'snrThreshold')
                    params.snrThreshold{i} = 0; % default: 0
                end

                if ~isfield(params,'useSortedSpikesKMeans')
                    [~,~,~,params.spikeChannelNumber{i},params.unitIDList{i}] = getMeanFiringRate(folderName,'electrodeList',params.electrodesFR{i},'goodPos',1,'timebin',params.FRtimebin,'snrThreshold',params.snrThreshold{i},'psthScreening','psthThreshold',params.psthThreshold{i},'useAllunitIDs','snrScreening','removeUnitID255',params.removeUnitID255{i},'useAsMultiunit',params.useAsMultiunit{i});
                else
                    if params.useSortedSpikesKMeans{i}
                        [~,~,~,params.spikeChannelNumber{i},params.unitIDList{i}] = getMeanFiringRate(folderName,'electrodeList',params.electrodesFR{i},'goodPos',1,'timebin',params.FRtimebin,'snrThreshold',params.snrThreshold{i},'psthScreening','psthThreshold',params.psthThreshold{i},'useAllunitIDs','useSortedSpikesKMeans','snrScreening','removeUnitID255',params.removeUnitID255{i},'useAsMultiunit',params.useAsMultiunit{i});
                    else
                        [~,~,~,params.spikeChannelNumber{i},params.unitIDList{i}] = getMeanFiringRate(folderName,'electrodeList',params.electrodesFR{i},'goodPos',1,'timebin',params.FRtimebin,'snrThreshold',params.snrThreshold{i},'psthScreening','psthThreshold',params.psthThreshold{i},'useAllunitIDs','snrScreening','removeUnitID255',params.removeUnitID255{i},'useAsMultiunit',params.useAsMultiunit{i});
                    end
                end

            end

            % Now generate the data for each condition
            for c = 1:numCols
                for r = 1:numRows
                    disp(['***************** processing: protocol ' num2str(i) '/' num2str(numProtocols) ', col ' num2str(c) '/' num2str(numCols) ', row ' num2str(r) '/' num2str(numRows)]);

                    % get the goodPos for this particular condition
                    if strncmpi(protocol,'GRF',3)
                        goodPos{r,c,i} = getGoodPosGRF(folderName,squeeze(valsMatrix(r,c,:)),removeCommonBadtrials);
                    elseif strncmpi(protocol,'CRS',3)
                        goodPos{r,c,i} = getGoodPosCRS(folderName,squeeze(valsMatrix(r,c,:)),removeCommonBadtrials);
                    end
                    if isfield(params,'ignoreStimPos')
                        goodPos{r,c,i} = setdiff(goodPos{r,c,i},params.ignoreStimPos{i});
                    end

                    % if badTrials are considered individually for each electrode
                    % then compute the goodPos for each electrode
                    numElecs = length(params.electrodesFR{i});
                    if removeIndividualBadtrials
                        for en = 1:numElecs
                            if strcmpi(params.reftype{1},'single')
                                params.elecbadtrialsFR{i,en} = allBadTrials{i}{params.electrodesFR{i}(en)};
                            else % bipolar - not worked out for LFP yet
                                params.elecbadtrialsFR{i,en} = union(allBadTrials{i}{params.chanElecs(en,1)},allBadTrials{i}{params.chanElecs(en,2)});
                            end
                            goodPos{r,c,i,en} = setdiff(goodPos{r,c,i},params.elecbadtrialsFR{i,en});
                            params.numRepeatsFR{r,c,i,en} = size(goodPos{r,c,i,en},2); % store the number of repeats for this condition
                        end
                    else
                        params.numRepeatsFR{r,c,i} = size(goodPos{r,c,i},2); % store the number of repeats for this condition
                    end

                    if params.fr
                        if ~isfield(params,'useSortedSpikesKMeans')
                            [spikeData{r,c,i},psthVals{r,c,i},params.FRxs,params.spikeChannelNumber{r,c,i},params.unitIDList{r,c,i}] = getMeanFiringRate(folderName,'electrodeList',params.spikeChannelNumber{i},'goodPos',squeeze(goodPos(r,c,i,:)),'timebin',params.FRtimebin,'unitIDList',params.unitIDList{i},'removeUnitID255',params.removeUnitID255{i},'useAsMultiunit',params.useAsMultiunit{i});
                        else
                            if params.useSortedSpikesKMeans{i}
                                [spikeData{r,c,i},psthVals{r,c,i},params.FRxs,params.spikeChannelNumber{r,c,i},params.unitIDList{r,c,i}] = getMeanFiringRate(folderName,'electrodeList',params.spikeChannelNumber{i},'goodPos',squeeze(goodPos(r,c,i,:)),'timebin',params.FRtimebin,'unitIDList',params.unitIDList{i},'useSortedSpikesKMeans','removeUnitID255',params.removeUnitID255{i},'useAsMultiunit',params.useAsMultiunit{i});
                            else
                                [spikeData{r,c,i},psthVals{r,c,i},params.FRxs,params.spikeChannelNumber{r,c,i},params.unitIDList{r,c,i}] = getMeanFiringRate(folderName,'electrodeList',params.spikeChannelNumber{i},'goodPos',squeeze(goodPos(r,c,i,:)),'timebin',params.FRtimebin,'unitIDList',params.unitIDList{i},'removeUnitID255',params.removeUnitID255{i},'useAsMultiunit',params.useAsMultiunit{i});
                            end
                        end
                    end
                end
            end
        end
    end
end

if ~savedDataExistsFR && params.saveDataFlag && processFRData % just load data if available and not generate the electrode FR data everytime
    disp('<><><><>~~ Saving goodPos FR data ~~<><><><> ');
    save(fullfile(saveDataFolder,tagFR),'spikeData','psthVals','params')
%     numRepeats = params.numRepeats;
%     save('elecData.mat','elecData','numRepeats');
end

params.valsMatrix = valsMatrix;
if isfield(params,'hMsgBox')
    message = params.message;
    hMsgBox = params.hMsgBox;
end
params.maxPlotsAlongX = 6;

if ~isfield(params,'drawplotsERP'); params.drawplotsERP = 0; end
if ~isfield(params,'drawplotsPSD'); params.drawplotsPSD = 0; end
if ~isfield(params,'drawplotsTF'); params.drawplotsTF = 0; end
if ~isfield(params,'drawplotsFR'); params.drawplotsFR = 0; end
if ~isfield(params,'newfiguresflag'); params.newfiguresflag = 0; end
if ~isfield(params,'returnERPData'); params.returnERPData = 1; end
if ~isfield(params,'returnPSDData'); params.returnPSDData = 1; end
if ~isfield(params,'returnTFData'); params.returnTFData = 1; end
if ~isfield(params,'returnFRData'); params.returnFRData = 1; end

% make colors
if strcmpi(params.colorstyle,'manual')
    params.colorsListType = 'singleSequence';
elseif strcmpi(params.colorstyle,'hsv')
    for r=1:numRows
        for c = 1:numCols
            params.colorsList{r,c} = hsv2rgb([valsMatrix{r,c,5}/360 valsMatrix{r,c,4} valsMatrix{r,c,6}/100]); 
            % make an hsv color value using the ori, sf and con values for 
            % this case.
            % [valsMatrix{r,c,5}/360 valsMatrix{r,c,4} valsMatrix{r,c,6}/100] => [ori sf con] => [h s v]
            if params.colorsList{r,c}==[0 0 0] % black
                params.colorsList{r,c} = [0.15 0.15 0.15]; % make it slightly greyish
            end
            if params.colorsList{r,c}==[1 1 1] % white
                params.colorsList{r,c} = [0.9 0.9 0.9]; % make it slightly greyish
            end
        end
    end
    params.colorsListType = 'matrix';
else
    vertices{1} = [1 0 0]; vertices{2} = [0 0 1]; vertices{3} = [0.5 1 0]; vertices{4} = [0.2 0.5 0.9]; vertices{5} = [0.3 0.4 0.1];
    boundaries = [0 0.4 0.7 0.95 1];
    if numRows*numCols==1
        params.colorsList = [0 0 0]; % only one color - black
    else
        params.colorsList = customColormap(numRows*numCols,vertices,boundaries);
    end
    params.colorsListType = 'singleSequence';
end

if params.showsubjectwise
    uniqueSubjects = unique(params.subjectslist); % find the list of unique subjects
    for us = 1:length(uniqueSubjects)
        subjectIndices(us,:) = find(ismember(params.subjectslist,uniqueSubjects)); % and the protocol indices corresponding to each subject
    end
else
    subjectIndices = 1:length(params.subjectslist); % otherwise, all subjects treated as one
end

for si = 1:size(subjectIndices,1) % sweep across subjects
    if params.erp
        newMessage = '>:>:>:>: plotting erps, erp measures :<:<:<:<';
        disp(newMessage);
        if isfield(params,'hMsgBox')
            message = updateMessageBox(hMsgBox,newMessage,message);
        end
        [~,dataOut.erpdata] = computeAndPlotMeanERPMeasures(elecData(:,:,subjectIndices(si,:)),params,params.drawplotsERP,params.newfiguresflag,params.returnERPData);
    end
    
    if params.psd
        newMessage = '>:>:>:>: plotting psd, psd measures :<:<:<:<';
        disp(newMessage);
        if isfield(params,'hMsgBox')
            message = updateMessageBox(hMsgBox,newMessage,message);
        end
        [~,dataOut.psddata] = computeAndPlotMeanPSDMeasures(elecData(:,:,subjectIndices(si,:)),params,params.drawplotsPSD,params.newfiguresflag,params.returnPSDData);
    end
    
    if params.tf
        newMessage = '>:>:>:>: plotting TF spectra, tf measures :<:<:<:<';
        disp(newMessage);
        if isfield(params,'hMsgBox')
            message = updateMessageBox(hMsgBox,newMessage,message);
        end
        if strcmpi(params.baselinestrategy,'common') && strcmpi(params.tfstrategy,'default')
            disp('*******common baseline, log at every electrode********');
            [~,dataOut.tfdata] = computeAndPlotMeanTFMeasuresDefault(elecData(:,:,subjectIndices(si,:)),params,params.drawplotsTF,params.newfiguresflag,params.returnTFData);
        else
            [~,dataOut.tfdata] = computeAndPlotMeanTFMeasures(elecData(:,:,subjectIndices(si,:)),params,params.drawplotsTF,params.newfiguresflag,params.returnTFData);
        end
    end
    
    if params.fr
        newMessage = '>:>:>:>: plotting Firing Rates, FR measures :<:<:<:<';
        disp(newMessage);
        if isfield(params,'hMsgBox')
            message = updateMessageBox(hMsgBox,newMessage,message);
        end
        [~,dataOut.frdata] = computeAndPlotMeanFRMeasures(spikeData(:,:,subjectIndices(si,:)),params,params.drawplotsFR,params.newfiguresflag,params.returnFRData,'psth',psthVals(:,:,subjectIndices(si,:)));
    end
    
    if params.stim
        newMessage = '>:>:>:>: drawing STIMULI :<:<:<:<';
        disp(newMessage);
        if isfield(params,'hMsgBox')
            message = updateMessageBox(hMsgBox,newMessage,message);
        end
        drawColorStimuliFromValsMatrix(params);
    end
end

end

function [allBadTrials, badTrials, badElecs, nameElec] = loadBadTrials(badTrialFile)
load(badTrialFile);
if ~exist('badElecs','var')
    badElecs = [];
end
end
