% This function is used to extract the digital data for the CRS Protocol.

% These are the following modes in which the CRSMap protocol has been used 
% so far.

% 1. The task is run in the fixation mode, in which target stimulus is off
% and all trials are catch trials. 
% Three gabors (gabor0/1/2 => S/R/C) are drawn and Digital codes are sent 
% for all parameters of S (i.e gabor0) and non-redundant parameters of
% gabor1 and gabor2 (R and C). Set hideSomeDigitalCode to 1 in this case.
% The target is assumed to be synchronous with the mapping stimulus 
% (gabor0) in this case.

% 2. The task is run in the fixation mode, in which target stimulus is off
% and all trials are catch trials. 
% gabor0 (i.e. S) is a null gabor and not drawn.
% Two gabors (gabor1/2 => R/C) are drawn and Digital codes are sent 
% for all non-redundant parameters of gabor1 and gabor2 (R and C). 
% Set hideSomeDigitalCode to 1 in this case.
% The target is assumed to be synchronous with the mapping stimulus 
% (gabor0) in this case.

% The mode depends on the specific protocol used (which is coded and sent
% in a variable Digital code - 'protocolNumber'. This can be read and used 
% to reconstruct the stimulus parameter values from the sent non-redundant
% Digital Codes)

function [goodStimNums,goodStimTimes,side] = extractDigitalDataCRS(folderExtract,ignoreTargetStimFlag,frameRate)

if ~exist('ignoreTargetStimFlag','var');   ignoreTargetStimFlag=0;      end
if ~exist('frameRate','var');              frameRate=100;               end

hideSomeDigitalCode=1; % Vinay - set this to 1 if you have sent partial code based on the protocol

stimResults = readDigitalCodesCRS(folderExtract,frameRate,hideSomeDigitalCode); % writes stimResults and trialResults
side = stimResults.side;
[goodStimNums,goodStimTimes] = getGoodStimNumsCRS(folderExtract,ignoreTargetStimFlag); % Good stimuli
save(fullfile(folderExtract,'goodStimNums.mat'),'goodStimNums');
end

% GRF Specific protocols
function [stimResults,trialResults,trialEvents] = readDigitalCodesCRS(folderOut,frameRate,hideSomeDigitalCode)

if ~exist('frameRate','var');              frameRate=100;               end
kForceQuit=7;

% Get the values of the following trial events for comparison with the dat
% file from lablib
trialEvents{1} = 'TS'; % Trial start
trialEvents{2} = 'TE'; % Trial End

load(fullfile(folderOut,'digitalEvents.mat'));

allDigitalCodesInDec = [digitalCodeInfo.codeNumber];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the times and values of the events in trialEvents

for i=1:length(trialEvents)
    pos = find(convertStrCodeToDec(trialEvents{i})==allDigitalCodesInDec);
    if isempty(pos)
        disp(['Code ' trialEvents{i} ' not found!!']);
    else
        trialResults(i).times = [digitalCodeInfo(pos).time]; %#ok<*AGROW>
        trialResults(i).value = [digitalCodeInfo(pos).value];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulus properties
azimuth          = [digitalCodeInfo(find(convertStrCodeToDec('AZ')==allDigitalCodesInDec)).value];
elevation        = [digitalCodeInfo(find(convertStrCodeToDec('EL')==allDigitalCodesInDec)).value];
contrast         = [digitalCodeInfo(find(convertStrCodeToDec('CO')==allDigitalCodesInDec)).value];
temporalFrequency= [digitalCodeInfo(find(convertStrCodeToDec('TF')==allDigitalCodesInDec)).value];
radius           = [digitalCodeInfo(find(convertStrCodeToDec('RA')==allDigitalCodesInDec)).value];
sigma            = [digitalCodeInfo(find(convertStrCodeToDec('SI')==allDigitalCodesInDec)).value];
spatialFrequency = [digitalCodeInfo(find(convertStrCodeToDec('SF')==allDigitalCodesInDec)).value];
orientation      = [digitalCodeInfo(find(convertStrCodeToDec('OR')==allDigitalCodesInDec)).value];
spatialPhase     = [digitalCodeInfo(find(convertStrCodeToDec('SP')==allDigitalCodesInDec)).value]; % [Vinay] - added spatial phase because it is being used in CRS

% Get timing
trialStartTimes = [digitalCodeInfo(find(convertStrCodeToDec('TS')==allDigitalCodesInDec)).time];
taskGaborTimes  = [digitalCodeInfo(find(convertStrCodeToDec('TG')==allDigitalCodesInDec)).time];
mapping0Times   = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time];
mapping1Times   = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time];
mapping2Times   = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]; % [Vinay] - added for CRS: the centre gabor
% [Vinay]: In CRS - S = gabor0, R = gabor1, C = gabor2
numTrials = length(trialStartTimes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In CRS, partial digital codes may be sent i.e. only those that are
% required to identify the stimulus completely, avoiding redundant codes.
% For eg. in some protocols the properties of Surround and Centre are
% matched (eg. same contrast, or orientation etc.). In such cases one may
% send values for only one of the two gabors. Additionally the
% 'protocolNumber' is sent as a digital code and this number identifies the
% specific protocol that is run. We use this information to fill all the
% parameter values for all the gabors.


% [Vinay] - Read the protocolNumber
stimResults.protocolNumber = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('PN')==allDigitalCodesInDec)).value]);

%----
% Adjust the stimResults.parameter based on the protocolNumber if partial
% code was sent.
% Any parameter is sent as a sequential train of values for gabor0, 1 and
% 2, for every stimulus. Therefore the values are in repeated groups of 1,
% 2 or 3 depending on the number of gabors that the values for that
% particular parameter are sent. For eg. if the values are sent for all
% three gabors they will be ordered as - [g0 g1 g2 g0 g1 g2 ...], if only
% the ring (g1) and centre (g2) values are sent then they would be - [g1 g2
% g1 g2 ...], if only surround values are sent it would be - [g0 g0 g0 ...]
% and so on. We take the values that are sent and repeat/arrange them
% so as to fill the values for all gabors before storing them.
% For eg. in Ring protocol all values for S (gabor0) are sent and only the
% 'radius' values are sent for R and C. Therefore radius values will be
% ordered as [g0 g1 g2 g0 g1 g2 ..] whereas all others are [g0 g0 g0 ..].
% Hence, radius remains unchanged, while for others the values sent for g0
% are repeated for g1 and g2 and the new parameter matrix obtained. For
% contrast, g0 values are repeated only for g2 whereas g1 is assigned a 
% value of 0 because Ring is always shown at 0% contrast in Ring Protocol.
% Similar procedure used for other protocols as per the specifics of each
% protocol
% NOTE: All digital codes are sent for protocolNumber 0 (noneProtocol), 9
% (cross-orientation Protocol). Therefore no re-arrangement is required in
% those cases
%------
if hideSomeDigitalCode
    switch stimResults.protocolNumber
        case 1 % Ring Protocol
            azi = [azimuth,azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation,elevation]; elevation = reshape(ele',[],1);
            ori = [orientation,orientation,orientation]; orientation = reshape(ori',[],1);
            tf = [temporalFrequency,temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            sp = [spatialPhase,spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
            con = [contrast,zeros(length(contrast),1),contrast]; contrast = reshape(con',[],1);
        case 2 % Contrast Ring Protocol
            azi = [azimuth,azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation,elevation]; elevation = reshape(ele',[],1);
            ori = [orientation,orientation,orientation]; orientation = reshape(ori',[],1);
            tf = [temporalFrequency,temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            sp = [spatialPhase,spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
            g1 = 1:2:length(contrast); g2 = 2:2:length(contrast);
            con = [contrast(g1),contrast(g2),contrast(g1)]; contrast = reshape(con',[],1);
        case 3 % Dual Contrast Protocol
            azi = [azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation]; elevation = reshape(ele',[],1);
            ori = [orientation,orientation]; orientation = reshape(ori',[],1);
            tf = [temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            sp = [spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
        case 4 % Dual Orientation Protocol
            azi = [azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation]; elevation = reshape(ele',[],1);
            con = [contrast,contrast]; contrast = reshape(con',[],1);
            tf = [temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            sp = [spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
        case 5 % Dual Phase Protocol
            azi = [azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation]; elevation = reshape(ele',[],1);
            con = [contrast,contrast]; contrast = reshape(con',[],1);
            tf = [temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            ori = [orientation,orientation]; orientation = reshape(ori',[],1);
        case 6 % Orientation Ring Protocol
            azi = [azimuth,azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation,elevation]; elevation = reshape(ele',[],1);
            con = [contrast,contrast,contrast]; contrast = reshape(con',[],1);
            tf = [temporalFrequency,temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            sp = [spatialPhase,spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
            g1 = 1:2:length(orientation); g2 = 2:2:length(orientation);
            ori = [orientation(g1),orientation(g2),orientation(g1)]; orientation = reshape(ori',[],1);
        case 7 % Phase Ring Protocol
            azi = [azimuth,azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation,elevation]; elevation = reshape(ele',[],1);
            con = [contrast,contrast,contrast]; contrast = reshape(con',[],1);
            tf = [temporalFrequency,temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            ori = [orientation,orientation,orientation]; orientation = reshape(ori',[],1);
            g1 = 1:2:length(spatialPhase); g2 = 2:2:length(spatialPhase);
            sp = [spatialPhase(g1),spatialPhase(g2),spatialPhase(g1)]; spatialPhase = reshape(sp',[],1);
        case 8 % Drifting Phase Protocol
            azi = [azimuth,azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation,elevation]; elevation = reshape(ele',[],1);
            con = [contrast,contrast,contrast]; contrast = reshape(con',[],1);
            sp = [spatialPhase,spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
            sig = [sigma,sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            ori = [orientation,orientation,orientation]; orientation = reshape(ori',[],1);
            g1 = 1:2:length(temporalFrequency); g2 = 2:2:length(temporalFrequency);
            tf = [temporalFrequency(g1),temporalFrequency(g2),temporalFrequency(g1)]; temporalFrequency = reshape(tf',[],1);
        case 10 % Annulus Fixed Protocol
            azi = [azimuth,azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation,elevation]; elevation = reshape(ele',[],1);
            ori = [orientation,orientation,orientation]; orientation = reshape(ori',[],1);
            tf = [temporalFrequency,temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            sp = [spatialPhase,spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
            g1 = 1:2:length(contrast); g2 = 2:2:length(contrast);
            con = [contrast(g1),contrast(g2),contrast(g1)]; contrast = reshape(con',[],1);
        case 11 % Colour Dual Protocol
            azi = [azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation]; elevation = reshape(ele',[],1);
            tf = [temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma]; sigma = reshape(sig',[],1);
            sp = [spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
        case 12 % Colour Ring Protocol
            azi = [azimuth,azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation,elevation]; elevation = reshape(ele',[],1);
            o1 = 1:2:length(orientation); o2 = 2:2:length(orientation); 
            ori = [orientation(o1),orientation(o2),orientation(o1)]; orientation = reshape(ori',[],1);
            tf = [temporalFrequency,temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma,sigma]; sigma = reshape(sig',[],1);
            sf1 = 1:2:length(spatialFrequency); sf2 = 2:2:length(spatialFrequency); 
            sf = [spatialFrequency(sf1),spatialFrequency(sf2),spatialFrequency(sf1)]; spatialFrequency = reshape(sf',[],1);
            sp = [spatialPhase,spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
            g1 = 1:2:length(contrast); g2 = 2:2:length(contrast);
            con = [contrast(g1),contrast(g2),contrast(g1)]; contrast = reshape(con',[],1);
    end
end


if isempty(azimuth) || isempty(elevation) || isempty(contrast) || isempty(temporalFrequency) ...
        || isempty(radius) || isempty(sigma) || isempty(spatialFrequency) || isempty(orientation) ...
        || isempty(spatialPhase)
    
    disp('Digital codes for the stimuli are not sent. Read from Lablib data file later ...');
else
    
    % Check the default case - only mapping0/1 is on, and only its stimulus properties are put out.
    
    if (max(diff([length(azimuth) length(elevation) length(contrast) length(temporalFrequency) ...
            length(radius) length(sigma) length(spatialFrequency) length(orientation) ...
            length(spatialPhase)])) > 0 )
        
        error('Length of stimulus properties are not even');
    else
        if((length(azimuth) == length(mapping0Times)) && isempty(mapping1Times) && isempty(mapping2Times))
            disp('Only Mapping 0 is used i.e. S ON, C & R OFF'); % Vinay, CRS: S ON, C & R OFF
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
            stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
            stimResults.side=0;

        elseif((length(azimuth) == length(mapping1Times)) && isempty(mapping0Times) && isempty(mapping2Times))
            disp('Only Mapping 1 is used i.e. R ON, C & S OFF'); % Vinay, CRS: R ON, C & S OFF
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
            stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
            stimResults.side=1;

        elseif((length(azimuth) == length(mapping2Times)) && isempty(mapping0Times) && isempty(mapping2Times))
            disp('Only Mapping 2 is used i.e. C ON, R & S OFF'); % Vinay, CRS: C ON, R & S OFF
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
            stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
            stimResults.side=2;

        elseif((length(azimuth) == (length(mapping0Times) + length(mapping1Times))) && isempty(mapping2Times))
            disp('Only Mapping 0 and 1 are used i.e S & R ON, C OFF'); % Vinay, CRS: S & R ON, C OFF 
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
            stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
            stimResults.side=[0 1];

        elseif((length(azimuth) == (length(mapping0Times) + length(mapping2Times))) && isempty(mapping1Times))
            disp('Only Mapping 0 and 2 are used i.e. S & C ON, R OFF'); % Vinay, CRS: S & C ON, R OFF
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
            stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
            stimResults.side=[0 2];

        elseif((length(azimuth) == (length(mapping1Times) + length(mapping2Times))) && isempty(mapping0Times))
            disp('Only Mapping 1 and 2 are used i.e. R & C ON, S OFF'); % Vinay, CRS: R & C ON, S OFF
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
            stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
            stimResults.side=[1 2]; 

        else
            disp('Digital codes from all gabors!!!');
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
            stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
            stimResults.side=[0 1 2];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instruction trials
instructionTrials = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('IT')==allDigitalCodesInDec)).value])';
if length(instructionTrials) ~= numTrials 
    error('Number of instruction trial entries different from numTrials');
end

% Catch trials
catchTrials = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('CT')==allDigitalCodesInDec)).value])';
if length(catchTrials) ~= numTrials 
    error('Number of catch trial entries different from numTrials');
end

% TrialCertify & TrialEnd (eotCode)
% These two entries may be repeated twice during force quit
trialCertify = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TC')==allDigitalCodesInDec)).value])';
eotCodes = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TE')==allDigitalCodesInDec)).value])';

forceQuits = find(eotCodes==kForceQuit);
numForceQuits = length(forceQuits);

% Vinay - it seems that when the task is stopped it gives an extra trialEnd
% code. Have to snub this, else it gives an error.
trialEndError = 0; % initialize
if ((length(trialResults(2).times) - length(trialResults(1).times)) == 1)
    trialResults(2).times(length(trialResults(2).times)) = [];
    trialResults(2).value(length(trialResults(2).times)) = [];
    trialEndError = 1;
end

% [Vinay] - if the final code of 'eotCodes' is 5 (False Alarm/quit) then it
% might be becuase of the stipulated blocks getting exhausted and the
% experiment stopping as a result of that. Check for this condition and if
% it is found to be so then truncate the last eotCode and trialCertify.
% Also, it seems that when the task is stopped it gives an extra trialEnd
% code if the blocks are not exhausted. So the last eotCodes and
% trialCertify have to be dropped. This has been taken care of in the
% trialResults reading code earlier and the trialEndError indicates if this
% error occured or not
if ((eotCodes(length(eotCodes))==5) || (trialEndError == 1))
    eotCodes(length(eotCodes)) = [];
    trialCertify(length(eotCodes)) = [];
end

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
numStimTask  = getStimPosPerTrial(trialStartTimes, taskGaborTimes);
numStimMap0  = getStimPosPerTrial(trialStartTimes, mapping0Times);
numStimMap1  = getStimPosPerTrial(trialStartTimes, mapping1Times);
numStimMap2  = getStimPosPerTrial(trialStartTimes, mapping2Times); % Vinay - for CRS

% [Vinay] - store the stim times of the C,R,S gabors
stimResults.time0 = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time]';
stimResults.time1 = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time]';
stimResults.time2 = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]';
% [Vinay] - assign numStims as any of numStimMap0/1/2 based on whichever
% is used. The numbers for the gabors which are drawn should be equal,
% so it doesn't matter which one is assigned to numStims. numStims says
% how many stimuli were shown on each of the trials

% Check if Task and Mapping stim nums are the same for non-intruction trials
nonInstructionTrials = find(instructionTrials==0);
if (max(abs(numStimTask(nonInstructionTrials) - numStimMap0(nonInstructionTrials)))==0)
    disp('Mapping0 and Task times are the same');
    numStimTask=numStimMap0;
    numStims = numStimMap0;
    stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time]';
    % stimResults.side = 0;
    taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TG')==allDigitalCodesInDec)).value])'; 
    
    if sum(taskType)==0 % Target is always null
        taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).value])';
    end
    
elseif (max(abs(numStimTask(nonInstructionTrials) - numStimMap1(nonInstructionTrials)))==0)
    disp('Mapping1 and Task times are the same');
    numStims = numStimMap1;
    stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time]';
    % stimResults.side = 1;
    taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TG')==allDigitalCodesInDec)).value])';
    
    if sum(taskType)==0 % Target is always null
        taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).value])';
    end

% [Vinay] - check gabor2, for CRS ----
elseif (max(abs(numStimTask(nonInstructionTrials) - numStimMap2(nonInstructionTrials)))==0)
    disp('Mapping2 and Task times are the same');
    numStims = numStimMap2;
    stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]';
    % stimResults.side = 1;
    taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TG')==allDigitalCodesInDec)).value])';
    
    if sum(taskType)==0 % Target is always null
        taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).value])';
    end
    
else
    disp('Mapping0/1/2 and Task times not the same');
    
    % [Vinay] - define another variable 'time' in stimResults to
    % store the onset time of stimulus display. This variable is
    % assigned the times corresponding to gabor0(surround) if it
    % is drawn (since it is the fist one to be drawn), else the
    % times corresponding to gabor1 (ring) and if even gabor1 is
    % undrawn then the times of gabor2
    
    if sum(numStimMap0)>0
        disp('Using Mapping0 times instead of task times...');
        numStimTask=numStimMap0;        % Assume task time is the same as mapping stimulus time                            
        numStims = numStimMap0;
        stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time]';
        
        % stimResults.side = 0;
        taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).value])'; % Assume task times and types are the same as M0
    
    elseif sum(numStimMap1)>0
        disp('Using Mapping1 times instead of task times...');
        numStimTask=numStimMap1;        % Assume task time is the same as mapping stimulus time
        numStims = numStimMap1;
        stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time]';
        
        % stimResults.side = 1;
        taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).value])'; % Assume task times and types are the same as M1
        
    elseif sum(numStimMap2)>0
        disp('Using Mapping2 times instead of task times...');
        numStimTask=numStimMap2;        % Assume task time is the same as mapping stimulus time
        numStims = numStimMap2;
        stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]';
        
        % stimResults.side = 2;
        taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).value])'; % Assume task times and types are the same as M2
    end
end

posTask = 0;
pos=0;
for i=1:numTrials
    taskTypeThisTrial = taskType(posTask+1:posTask+numStims(i)); 
    % [Vinay] - numStimTask works here only if there are no instruction 
    % trials, otherwise numStims contains the number of stimuli in each 
    % trial and works for all cases
    if (numStims(i)>0)
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
