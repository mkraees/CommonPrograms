% This function is used to extract the digital data for the GRF Protocol.

% These are the following modes in which the GaborRFMap protocol has been used so far.

% 1. Target and mapping stimulus 0 are on and are presented synchronously,
% while mapping stimulus 1 is off. Digital codes are sent only for map0 stimulus.

% 2. The task is run in the fixation mode, in which target stimulus is off
% and all trials are catch trials. Digital codes are sent only for map0
% stimulus. The invisible target is assumed to be synchronous with the
% mapping stimulus in this case. 

% This is copied from extractDigitalDataGRF. This reads all the digital
% data from LL file

function [goodStimNums,goodStimTimes,side] = extractDigitalDataGRFLL(folderExtract,ignoreTargetStimFlag,frameRate)

if ~exist('ignoreTargetStimFlag','var');   ignoreTargetStimFlag=0;      end
if ~exist('frameRate','var');              frameRate=100;               end

stimResults = readDigitalCodesGRF(folderExtract,frameRate); % writes stimResults and trialResults
side = stimResults.side;
[goodStimNums,goodStimTimes] = getGoodStimNumsGRF(folderExtract,ignoreTargetStimFlag); % Good stimuli
save(fullfile(folderExtract,'goodStimNums.mat'),'goodStimNums','goodStimTimes');
end

% GRF Specific protocols
function [stimResults,trialResults,trialEvents] = readDigitalCodesGRF(folderExtract,frameRate)

if ~exist('frameRate','var');              frameRate=100;               end
kForceQuit=7;

% TrialEvents are actually not useful - just keeping for compatibility with
% older programs.

trialEvents{1} = 'TS'; % Trial start
trialEvents{2} = 'TE'; % Trial End

load(fullfile(folderExtract,'digitalEvents.mat'));

allDigitalCodesInDec = [digitalCodeInfo.codeNumber];
useSingelITC18Flag=1;
if max(allDigitalCodesInDec)<=128
    useSimpleCodeFlag=1;
else
    useSimpleCodeFlag=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the times and values of the events in trialEvents

for i=1:length(trialEvents)
    pos = find(convertStrCodeToDec(trialEvents{i},useSingelITC18Flag,useSimpleCodeFlag)==allDigitalCodesInDec);
    if isempty(pos)
        disp(['Code ' trialEvents{i} ' not found!!']);
    else
        trialResults(i).times = [digitalCodeInfo(pos).time]; %#ok<*AGROW>
        trialResults(i).value = [digitalCodeInfo(pos).value];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Lablib data structure
load(fullfile(folderExtract,'LL.mat'));

if sum(LL.stimType1)>0
    activeSide=0;
elseif sum(LL.stimType2)>0
    activeSide=1;
end

%%%%%%%%%%%%%%%%% Get info from LL to construct stimResults %%%%%%%%%%%%%%%
if activeSide==0 % Map0
    validMap = find(LL.stimType1>0);
    aziLL = LL.azimuthDeg1(validMap);
    eleLL = LL.elevationDeg1(validMap);
    sigmaLL = LL.sigmaDeg1(validMap);
    
    if isfield(LL,'radiusDeg1')
        radiusExists = 1;
        radiusLL = LL.radiusDeg1(validMap);
    else
        radiusExists = 0;
    end
    sfLL = LL.spatialFreqCPD1(validMap);
    oriLL = LL.orientationDeg1(validMap);
    conLL = LL.contrastPC1(validMap); 
    tfLL = LL.temporalFreqHz1(validMap); 
    timeLL = LL.time1(validMap)/1000;
    mapping0Times = timeLL;
    taskType = LL.stimType1(validMap);
    
elseif activeSide==1 % Map2
    
    validMap = find(LL.stimType2>0);
    aziLL = LL.azimuthDeg2(validMap);
    eleLL = LL.elevationDeg2(validMap);
    sigmaLL = LL.sigmaDeg2(validMap);
    if isfield(LL,'radiusDeg2')
        radiusExists = 1;
        radiusLL = LL.radiusDeg2(validMap);
    else
        radiusExists = 0;
    end
    sfLL = LL.spatialFreqCPD2(validMap);
    oriLL = LL.orientationDeg2(validMap);
    conLL = LL.contrastPC2(validMap); 
    tfLL = LL.temporalFreqHz2(validMap);
    timeLL = LL.time2(validMap)/1000;
    mapping1Times = timeLL;
    taskType = LL.stimType2(validMap);
end

stimResults.azimuth = aziLL;
stimResults.elevation = eleLL;
stimResults.contrast = conLL;
stimResults.temporalFrequency = tfLL;
if radiusExists
    stimResults.radius = radiusLL;
end
stimResults.sigma = sigmaLL;
stimResults.orientation = oriLL;
stimResults.spatialFrequency = sfLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get timing
    
trialStartTimes = [digitalCodeInfo(find(convertStrCodeToDec('TS',useSingelITC18Flag,useSimpleCodeFlag)==allDigitalCodesInDec)).time];
eotCodes = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TE',useSingelITC18Flag,useSimpleCodeFlag)==allDigitalCodesInDec)).value])';
trialStartTimesLL = LL.startTime;
eotCodesLL = LL.eotCode;
if useSimpleCodeFlag
    eotCodes=eotCodesLL;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compare TS and TE data %%%%%%%%%%%%%%%%%%%%%%%
diffTD = diff(trialStartTimes); diffTL = diff(trialStartTimesLL);

maxDiffMS = 1000*max(abs(diffTD(:) - diffTL(:)));
dEOT = max(abs(diff(eotCodes(:)-eotCodesLL(:))));

maxDiffCutoffMS = 5; % throw an error if the difference exceeds 5 ms
if maxDiffMS > maxDiffCutoffMS || dEOT > 0
    error('The digital codes do not match with the LL data...');
else
    disp(['Maximum difference between LL and LFP/EEG start times: ' num2str(maxDiffMS) ' ms']);
end
numTrials = length(trialStartTimes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instruction trials
instructionTrials = LL.instructTrial;
if length(instructionTrials) ~= numTrials 
    error('Number of instruction trial entries different from numTrials');
end

% Catch trials
catchTrials = LL.catchTrial;
if length(catchTrials) ~= numTrials 
    error('Number of catch trial entries different from numTrials');
end

% TrialCertify & TrialEnd (eotCode)
% These two entries may be repeated twice during force quit
trialCertify = LL.trialCertify;

forceQuits = find(eotCodes==kForceQuit);
numForceQuits = length(forceQuits);

if length(eotCodes)-numForceQuits == numTrials
    disp(['numTrials: ' num2str(numTrials) ' numEotCodes: '  ...
        num2str(length(eotCodes)) ', ForceQuits: ' num2str(numForceQuits)]);
    goodEOTPos = find(eotCodes ~=kForceQuit);
    eotCodes = eotCodes(goodEOTPos);
    trialCertify = trialCertify(goodEOTPos);
else
     disp(['numTrials: ' num2str(numTrials) ' numEotCodes: '  ...
        num2str(length(eotCodes)) ', forcequits: ' num2str(numForceQuits)]);
    error('ForceQuit pressed after trial started'); % TODO - deal with this case
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numStims  = getStimPosPerTrial(trialStartTimesLL,timeLL);
stimResults.side = activeSide;

posTask = 0;
pos=0;
for i=1:numTrials 
    
    taskTypeThisTrial = taskType(posTask+1:posTask+numStims(i));
    if (numStims(i)>0)
        stimTimesFromTrialStart = timeLL(pos+1:pos+numStims(i)) - trialStartTimesLL(i);
        stimResults.time(pos+1:pos+numStims(i)) = trialStartTimes(i) + stimTimesFromTrialStart; % times relative to the LFP/EEG data, not LL Data
    
        stimResults.type(pos+1:pos+numStims(i)) = taskTypeThisTrial;
        stimResults.trialNumber(pos+1:pos+numStims(i)) = i;
        stimResults.stimPosition(pos+1:pos+numStims(i)) = 1:numStims(i);
        
        if stimResults.side==0
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
                (mapping0Times(pos+1:pos+numStims(i)) - mapping0Times(pos+1))*frameRate;
        elseif stimResults.side==1
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
                (mapping1Times(pos+1:pos+numStims(i)) - mapping1Times(pos+1))*frameRate;
        end
        
        stimResults.instructionTrials(pos+1:pos+numStims(i)) = instructionTrials(i); %always zero
        stimResults.catch(pos+1:pos+numStims(i)) = catchTrials(i);
        stimResults.eotCodes(pos+1:pos+numStims(i)) = eotCodes(i);
        stimResults.trialCertify(pos+1:pos+numStims(i)) = trialCertify(i);
        pos = pos+numStims(i);
    end
    posTask = posTask+numStims(i);
end

% Save in folderExtract
save(fullfile(folderExtract,'stimResults.mat'),'stimResults');
save(fullfile(folderExtract,'trialResults.mat'),'trialEvents','trialResults');

end
function [goodStimNums,goodStimTimes] = getGoodStimNumsGRF(folderExtract,ignoreTargetStimFlag)

if ~exist('ignoreTargetStimFlag','var');      ignoreTargetStimFlag=0;   end

load(fullfile(folderExtract,'stimResults.mat'));

totalStims = length(stimResults.eotCodes);
disp(['Number of trials: ' num2str(max(stimResults.trialNumber))]);
disp(['Number of stimuli: ' num2str(totalStims)]);

% exclude uncertified trials, catch trials and instruction trials
tc = find(stimResults.trialCertify==1);
it = find(stimResults.instructionTrials==1);
ct = find(stimResults.catch==1);

if ignoreTargetStimFlag
    badStimNums = [it tc]; % catch trials are now considered good
else
    badStimNums = [it tc ct];
end

%eottypes
% 0 - correct, 1 - wrong, 2-failed, 3-broke, 4-ignored, 5-False
% Alarm/quit, 6 - distracted, 7 - force quit
%disp('Analysing correct, wrong and failed trials');
%badEOTs = find(stimResults.eotCodes>2); 
%disp('Analysing correct and wrong trials')
%badEOTs = find(stimResults.eotCodes>1); 
disp('Analysing only correct trials')
badEOTs = find(stimResults.eotCodes>0); 
badStimNums = [badStimNums badEOTs];

goodStimNums = setdiff(1:totalStims,unique(badStimNums));

% stim types
% 0 - Null, 1 - valid, 2 - target, 3 - frontpadding, 4 - backpadding
if ~ignoreTargetStimFlag
    disp('Only taking valid stims ');
    validStims = find(stimResults.type==1);
    goodStimNums = intersect(goodStimNums,validStims);
    
    %%%%%%%%%%%%%% Remove bad stimuli after target %%%%%%%%%%%%%%%%%%%%
    
    clear trialNums stimPos
    trialNums = stimResults.trialNumber(goodStimNums);
    stimPos   = stimResults.stimPosition(goodStimNums);
    
    % Get the target positions of the trialNums
    clear goodTrials
    goodTrials = unique(trialNums);
    
    clear targetPos
    for i=1:length(goodTrials)
        allStimWithThisTrialNum = find(stimResults.trialNumber==goodTrials(i));
        
        if sum(stimResults.catch(allStimWithThisTrialNum))>0        % catch trials
            targetPos(trialNums==goodTrials(i)) = inf; %#ok<*AGROW>
        else
            targetPos(trialNums==goodTrials(i)) = find(stimResults.type(allStimWithThisTrialNum)==2);
        end
    end
    
    validStimuliAfterTarget = find(stimPos>targetPos);
    if ~isempty(validStimuliAfterTarget)
        disp([num2str(length(validStimuliAfterTarget)) ' out of ' num2str(length(goodStimNums)) ' stimuli after target']);
        save(fullfile(folderExtract,'validStimAfterTarget.mat'),'validStimuliAfterTarget');
    end
    
    goodStimNums(validStimuliAfterTarget)=[];
end
disp(['Number of good stimuli: ' num2str(length(goodStimNums))]);
goodStimTimes = stimResults.time(goodStimNums);
end
function outNum = convertUnits(num,f)

if ~exist('f','var');                       f=1;                        end

for i=1:length(num)
    if num(i) > 16384
        num(i)=num(i)-32768;
    end
end
outNum = num/f;
end
function [numStim,stimOnPos] = getStimPosPerTrial(trialStartTimes, stimStartTimes)

numTrials = length(trialStartTimes);

stimOnPos = cell(1,numTrials);
numStim   = zeros(1,numTrials);

for i=1:numTrials-1
    stimOnPos{i} = intersect(find(stimStartTimes>=trialStartTimes(i)),find(stimStartTimes<trialStartTimes(i+1)));
    numStim(i) = length(stimOnPos{i});
end
stimOnPos{numTrials} = find(stimStartTimes>=trialStartTimes(numTrials));
numStim(numTrials) = length(stimOnPos{numTrials});
end
