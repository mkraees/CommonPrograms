% This function is used to extract the digital data for the CRS Protocol.

% These are the following modes in which the CRSMap protocol has been used 
% so far.

% 1. The task is run in the fixation mode, in which target stimulus is off
% and all trials are catch trials. 
% Three gabors (gabor0/1/2 => S/R/C) are drawn
% The target is assumed to be synchronous with the mapping stimulus 
% (gabor0) in this case.

% 2. The task is run in the fixation mode, in which target stimulus is off
% and all trials are catch trials. 
% gabor0 (i.e. S) is a null gabor and not drawn.
% Two gabors (gabor1/2 => R/C) are drawn 
% The target is assumed to be synchronous with the mapping stimulus 
% (gabor0) in this case.

% The mode depends on the specific protocol used

function [goodStimNums,goodStimTimes,side] = extractDigitalDataCRSLL(folderExtract,ignoreTargetStimFlag,frameRate)

if ~exist('ignoreTargetStimFlag','var');   ignoreTargetStimFlag=0;      end
if ~exist('frameRate','var');              frameRate=100;               end

% hideSomeDigitalCode=1; % Vinay - set this to 1 if you have sent partial code based on the protocol
% hideSomeDigitalCode not required while extracting from LL data file

stimResults = readDigitalCodesCRS(folderExtract,frameRate); % writes stimResults and trialResults
side = stimResults.side;
[goodStimNums,goodStimTimes] = getGoodStimNumsCRS(folderExtract,ignoreTargetStimFlag); % Good stimuli
save(fullfile(folderExtract,'goodStimNums.mat'),'goodStimNums');
end

% GRF Specific protocols
function [stimResults,trialResults,trialEvents] = readDigitalCodesCRS(folderOut,frameRate)

if ~exist('frameRate','var');              frameRate=100;               end
kForceQuit=7;

% Get the values of the following trial events for comparison with the dat
% file from lablib
trialEvents{1} = 'TS'; % Trial start
trialEvents{2} = 'TE'; % Trial End

load(fullfile(folderOut,'digitalEvents.mat'));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Lablib data structure
load(fullfile(folderOut,'LL.mat'));

% Determine active gabors
s1 = sum(LL.stimType1); % This is 0 if the gabor is null
s2 = sum(LL.stimType2);
s3 = sum(LL.stimType3);
if (s1>0 && s2==0 && s3==0)
    activeSide=0;
elseif (s1==0 && s2>0 && s3==0)
    activeSide=1;
elseif (s1==0 && s2==0 && s3>0)
    activeSide=2;
elseif (s1>0 && s2>0 && s3==0)
    activeSide=[0 1];
elseif (s1>0 && s2==0 && s3>0)
    activeSide=[0 2];
elseif (s1==0 && s2>0 && s3>0)
    activeSide=[1 2];
else
    activeSide=[0 1 2];
end
stimResults.side = activeSide;

% [Vinay] - In CRS index 0 in LL represents the task gabor, whereas index
% 1, 2 and 3 represent gabor0 (S), gabor1 (R) and gabor2 (C) respectively
% activeSide cane be 0,1,2,[0 1],[0 2],[1 2] or [0 1 2]
% Irrespective of that, the parameter values for each gabor are stored.
% Therefore we just get the values for each one of them.
    
    validMap = find(LL.stimType1>0);
    aziLL1 = LL.azimuthDeg1(validMap);
    eleLL1 = LL.elevationDeg1(validMap);
    sigmaLL1 = LL.sigmaDeg1(validMap);
    radiusLL1 = LL.radiusDeg1(validMap);
    sfLL1 = LL.spatialFreqCPD1(validMap);
    oriLL1 = LL.orientationDeg1(validMap);
    conLL1 = LL.contrastPC1(validMap); 
    tfLL1 = LL.temporalFreqHz1(validMap);
    spLL1 = LL.SpatialPhaseDeg1(validMap);
    timeLL1 = LL.time1(validMap)/1000;
    mapping0Times = timeLL1;
    
    validMap = find(LL.stimType2>0);
    aziLL2 = LL.azimuthDeg2(validMap);
    eleLL2 = LL.elevationDeg2(validMap);
    sigmaLL2 = LL.sigmaDeg2(validMap);
    radiusLL2 = LL.radiusDeg2(validMap);
    sfLL2 = LL.spatialFreqCPD2(validMap);
    oriLL2 = LL.orientationDeg2(validMap);
    conLL2 = LL.contrastPC2(validMap); 
    tfLL2 = LL.temporalFreqHz2(validMap);
    spLL2 = LL.SpatialPhaseDeg2(validMap);
    timeLL2 = LL.time2(validMap)/1000;
    mapping1Times = timeLL2;
    
    validMap = find(LL.stimType3>0);
    aziLL3 = LL.azimuthDeg3(validMap);
    eleLL3 = LL.elevationDeg3(validMap);
    sigmaLL3 = LL.sigmaDeg3(validMap);
    radiusLL3 = LL.radiusDeg3(validMap);
    sfLL3 = LL.spatialFreqCPD3(validMap);
    oriLL3 = LL.orientationDeg3(validMap);
    conLL3 = LL.contrastPC3(validMap); 
    tfLL3 = LL.temporalFreqHz3(validMap);
    spLL3 = LL.SpatialPhaseDeg3(validMap);
    timeLL3 = LL.time3(validMap)/1000;
    mapping2Times = timeLL3;
    
    
    % Combine all values
    aziLL = [aziLL1; aziLL2; aziLL3]; aziLL = (reshape(aziLL,[],1))';
    eleLL = [eleLL1; eleLL2; eleLL3]; eleLL = (reshape(eleLL,[],1))';
    sigmaLL = [sigmaLL1; sigmaLL2; sigmaLL3]; sigmaLL = (reshape(sigmaLL,[],1))';
    sfLL = [sfLL1; sfLL2; sfLL3]; sfLL = (reshape(sfLL,[],1))';
    oriLL = [oriLL1; oriLL2; oriLL3]; oriLL = (reshape(oriLL,[],1))';
    conLL = [conLL1; conLL2; conLL3]; conLL = (reshape(conLL,[],1))';
    tfLL = [tfLL1; tfLL2; tfLL3]; tfLL = (reshape(tfLL,[],1))';
    spLL = [spLL1; spLL2; spLL3]; spLL = (reshape(spLL,[],1))';
    radiusLL = [radiusLL1; radiusLL2; radiusLL3]; radiusLL = (reshape(radiusLL,[],1))';
    

    stimResults.azimuth = aziLL;
    stimResults.elevation = eleLL;
    stimResults.contrast = conLL;
    stimResults.temporalFrequency = tfLL;
    stimResults.radius = radiusLL;
    stimResults.sigma = sigmaLL;
    stimResults.orientation = oriLL;
    stimResults.spatialFrequency = sfLL;
    stimResults.spatialPhase = spLL;

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
% [Vinay] - number of stim drawn with a gabor
numStimMap0  = getStimPosPerTrial(trialStartTimesLL, timeLL1); % S
numStimMap1  = getStimPosPerTrial(trialStartTimesLL, timeLL2); % R
numStimMap2  = getStimPosPerTrial(trialStartTimesLL, timeLL3); % C

% [Vinay] - store the stim times of the gabors
stimResults.time0 = timeLL1; % S
stimResults.time1 = timeLL2; % R
stimResults.time2 = timeLL3; % C

% [Vinay] - define another variable 'time' in stimResults to
% store the onset time of stimulus display. This variable is
% assigned the times corresponding to gabor0(surround) if it
% is drawn (since it is the fist one to be drawn), else the
% times corresponding to gabor1 (ring) and if even gabor1 is
% undrawn then the times of gabor2. Accordingly the taskType and numStims
% are also assigned in this order/preference. The numbers for the gabors 
% which are drawn should be equal, so it doesn't matter which one is 
% assigned to numStims. numStims says how many stimuli were shown on each 
% of the trials

if min(activeSide) == 0 % gabor0 i.e. S is drawn
    validMap = LL.stimType1>0;
    taskType = LL.stimType1(validMap);
    numStims = numStimMap0;
    timeLL = timeLL1;
elseif min(activeSide) == 1 % no gabor0, but gabor1 i.e. R is drawn
    validMap = LL.stimType2>0;
    taskType = LL.stimType2(validMap);
    numStims = numStimMap1;
    timeLL = timeLL2;
else % only gabor2 i.e. C is drawn
    validMap = LL.stimType3>0;
    taskType = LL.stimType3(validMap);
    numStims = numStimMap2;
    timeLL = timeLL3;
end


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
        elseif stimResults.side==2
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
                (mapping2Times(pos+1:pos+numStims(i)) - mapping2Times(pos+1))*frameRate;
        elseif min(stimResults.side)==1 % [1 2] case
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
                (mapping1Times(pos+1:pos+numStims(i)) - mapping1Times(pos+1))*frameRate;
        else % for all other cases gabor0 is drawn, so take its times i.e. for [0 1], [0 2], [0 1 2]
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
                (mapping0Times(pos+1:pos+numStims(i)) - mapping0Times(pos+1))*frameRate;
        end
        
        stimResults.instructionTrials(pos+1:pos+numStims(i)) = instructionTrials(i); %always zero
        stimResults.catch(pos+1:pos+numStims(i)) = catchTrials(i);
        stimResults.eotCodes(pos+1:pos+numStims(i)) = eotCodes(i);
        stimResults.trialCertify(pos+1:pos+numStims(i)) = trialCertify(i);
        pos = pos+numStims(i);
    end
    posTask = posTask+numStims(i);
end

% Save in folderOut
save(fullfile(folderOut,'stimResults.mat'),'stimResults');
save(fullfile(folderOut,'trialResults.mat'),'trialEvents','trialResults');

end
function [goodStimNums,goodStimTimes] = getGoodStimNumsCRS(folderOut,ignoreTargetStimFlag)

if ~exist('ignoreTargetStimFlag','var');      ignoreTargetStimFlag=0;   end

load(fullfile(folderOut,'stimResults.mat'));

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
        save(fullfile(folderOut,'validStimAfterTarget.mat'),'validStimuliAfterTarget');
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
