% This program generates the parameterCombinations variable from the
% stimResults
function parameterCombinations = getDisplayCombinationsGRF(folderOut,goodStimNums,radiusToSigmaRatio)

if ~exist('radiusToSigmaRatio','var');      radiusToSigmaRatio=3;       end

load(fullfile(folderOut,'stimResults.mat'));

% Five parameters are chosen:
% 1. Azimuth
% 2. Elevation
% 3. Sigma, Radius 
% 4. Spatial Frequency
% 5. Orientation
% 6. Contrast
% 7. Temporal Frequency

% Parameters index
parameters{1} = 'azimuth';
parameters{2} = 'elevation';
parameters{3} = 'sigma';
parameters{4} = 'spatialFrequency';
parameters{5} = 'orientation';
parameters{6} = 'contrast';
parameters{7} = 'temporalFrequency'; %#ok<NASGU>
parameters{8} = 'gaborNumber';

% get Contrast
aValsAll  = stimResults.azimuth;
eValsAll  = stimResults.elevation;
if isfield(stimResults,'radius')
    sValsAll  = stimResults.radius/radiusToSigmaRatio;
else
    sValsAll  = stimResults.sigma;
end
fValsAll  = stimResults.spatialFrequency;
oValsAll  = stimResults.orientation;
cValsAll  = stimResults.contrast;
tValsAll  = stimResults.temporalFrequency;

% Vinay - save parametersCombination with an extra dimension corresponding
% to the gabor number in the case when both gabors are drawn. We,
% therefore, check for the number of gabors and save accordingly.
numOfGabors = length(stimResults.side);

if numOfGabors ~= 1 % i.e when both left and right gabors are drawn
        
        aValsAllGabor = aValsAll;
        eValsAllGabor = eValsAll;
        sValsAllGabor = sValsAll;
        fValsAllGabor = fValsAll;
        oValsAllGabor = oValsAll;
        cValsAllGabor = cValsAll;
        tValsAllGabor = tValsAll;
    
    for gaborNum = 1:2
        
        gaborIndices = gaborNum:2:length(aValsAllGabor); % [Vinay] - index values corresponding to 'gaborNum' gabor
        
        aValsAll = aValsAllGabor(gaborIndices);
        eValsAll = eValsAllGabor(gaborIndices);
        sValsAll = sValsAllGabor(gaborIndices);
        fValsAll = fValsAllGabor(gaborIndices);
        oValsAll = oValsAllGabor(gaborIndices);
        cValsAll = cValsAllGabor(gaborIndices);
        tValsAll = tValsAllGabor(gaborIndices);
        
        aValsGood = aValsAll(goodStimNums);
        eValsGood = eValsAll(goodStimNums);
        sValsGood = sValsAll(goodStimNums);
        fValsGood = fValsAll(goodStimNums);
        oValsGood = oValsAll(goodStimNums);
        cValsGood = cValsAll(goodStimNums);
        tValsGood = tValsAll(goodStimNums);

        aValsUnique = unique(aValsGood); aLen = length(aValsUnique);
        eValsUnique = unique(eValsGood); eLen = length(eValsUnique);
        sValsUnique = unique(sValsGood); sLen = length(sValsUnique);
        fValsUnique = unique(fValsGood); fLen = length(fValsUnique);
        oValsUnique = unique(oValsGood); oLen = length(oValsUnique);
        cValsUnique = unique(cValsGood); cLen = length(cValsUnique);
        tValsUnique = unique(tValsGood); tLen = length(tValsUnique);
        
        % display
        disp(['Gabor ' num2str(gaborNum - 1)]); % Vinay - display kGabor number
        disp(['Number of unique azimuths: ' num2str(aLen)]);
        disp(['Number of unique elevations: ' num2str(eLen)]);
        disp(['Number of unique sigmas: ' num2str(sLen)]);
        disp(['Number of unique Spatial freqs: ' num2str(fLen)]);
        disp(['Number of unique orientations: ' num2str(oLen)]);
        disp(['Number of unique contrasts: ' num2str(cLen)]);
        disp(['Number of unique temporal freqs: ' num2str(tLen)]);

        % If more than one value, make another entry with all values
        if (aLen > 1);           aLen=aLen+1;                    end
        if (eLen > 1);           eLen=eLen+1;                    end
        if (sLen > 1);           sLen=sLen+1;                    end
        if (fLen > 1);           fLen=fLen+1;                    end
        if (oLen > 1);           oLen=oLen+1;                    end
        if (cLen > 1);           cLen=cLen+1;                    end
        if (tLen > 1);           tLen=tLen+1;                    end
        
        allPos = 1:length(goodStimNums);
        disp(['total combinations: ' num2str((aLen)*(eLen)*(sLen)*(fLen)*(oLen)*(cLen)*(tLen))]);
        
        for a=1:aLen
            if a==aLen
                aPos = allPos;
            else
                aPos = find(aValsGood == aValsUnique(a));
            end

            for e=1:eLen
                if e==eLen
                    ePos = allPos;
                else
                    ePos = find(eValsGood == eValsUnique(e));
                end

                for s=1:sLen
                    if s==sLen
                        sPos = allPos;
                    else
                        sPos = find(sValsGood == sValsUnique(s));
                    end

                    for f=1:fLen
                        if f==fLen
                            fPos = allPos;
                        else
                            fPos = find(fValsGood == fValsUnique(f));
                        end

                        for o=1:oLen
                            if o==oLen
                                oPos = allPos;
                            else
                                oPos = find(oValsGood == oValsUnique(o));
                            end

                            for c=1:cLen
                                if c==cLen
                                    cPos = allPos;
                                else
                                    cPos = find(cValsGood == cValsUnique(c));
                                end

                                for t=1:tLen
                                    if t==tLen
                                        tPos = allPos;
                                    else
                                        tPos = find(tValsGood == tValsUnique(t));
                                    end


                                    aePos = intersect(aPos,ePos);
                                    aesPos = intersect(aePos,sPos);
                                    aesfPos = intersect(aesPos,fPos);
                                    aesfoPos = intersect(aesfPos,oPos);
                                    aesfocPos = intersect(aesfoPos,cPos);
                                    aesfoctPos = intersect(aesfocPos,tPos);
                                    parameterCombinations{a,e,s,f,o,c,t,gaborNum} = aesfoctPos; %#ok<AGROW>
                                end
                            end
                        end
                    end
                end
            end
        end
        
    end
    save(fullfile(folderOut,'parameterCombinations.mat'),'parameters','parameterCombinations', ...
        'aValsUnique','eValsUnique','sValsUnique','fValsUnique','oValsUnique','cValsUnique','tValsUnique');
        
else
if ~isempty(aValsAll)

    aValsGood = aValsAll(goodStimNums);
    eValsGood = eValsAll(goodStimNums);
    sValsGood = sValsAll(goodStimNums);
    fValsGood = fValsAll(goodStimNums);
    oValsGood = oValsAll(goodStimNums);
    cValsGood = cValsAll(goodStimNums);
    tValsGood = tValsAll(goodStimNums);

    aValsUnique = unique(aValsGood); aLen = length(aValsUnique);
    eValsUnique = unique(eValsGood); eLen = length(eValsUnique);
    sValsUnique = unique(sValsGood); sLen = length(sValsUnique);
    fValsUnique = unique(fValsGood); fLen = length(fValsUnique);
    oValsUnique = unique(oValsGood); oLen = length(oValsUnique);
    cValsUnique = unique(cValsGood); cLen = length(cValsUnique);
    tValsUnique = unique(tValsGood); tLen = length(tValsUnique);

    % display
    disp(['Number of unique azimuths: ' num2str(aLen)]);
    disp(['Number of unique elevations: ' num2str(eLen)]);
    disp(['Number of unique sigmas: ' num2str(sLen)]);
    disp(['Number of unique Spatial freqs: ' num2str(fLen)]);
    disp(['Number of unique orientations: ' num2str(oLen)]);
    disp(['Number of unique contrasts: ' num2str(cLen)]);
    disp(['Number of unique temporal freqs: ' num2str(tLen)]);

    % If more than one value, make another entry with all values
    if (aLen > 1);           aLen=aLen+1;                    end
    if (eLen > 1);           eLen=eLen+1;                    end
    if (sLen > 1);           sLen=sLen+1;                    end
    if (fLen > 1);           fLen=fLen+1;                    end
    if (oLen > 1);           oLen=oLen+1;                    end
    if (cLen > 1);           cLen=cLen+1;                    end
    if (tLen > 1);           tLen=tLen+1;                    end

    allPos = 1:length(goodStimNums);
    disp(['total combinations: ' num2str((aLen)*(eLen)*(sLen)*(fLen)*(oLen)*(cLen)*(tLen))]);

    for a=1:aLen
        if a==aLen
            aPos = allPos;
        else
            aPos = find(aValsGood == aValsUnique(a));
        end

        for e=1:eLen
            if e==eLen
                ePos = allPos;
            else
                ePos = find(eValsGood == eValsUnique(e));
            end

            for s=1:sLen
                if s==sLen
                    sPos = allPos;
                else
                    sPos = find(sValsGood == sValsUnique(s));
                end

                for f=1:fLen
                    if f==fLen
                        fPos = allPos;
                    else
                        fPos = find(fValsGood == fValsUnique(f));
                    end

                    for o=1:oLen
                        if o==oLen
                            oPos = allPos;
                        else
                            oPos = find(oValsGood == oValsUnique(o));
                        end
                        
                        for c=1:cLen
                            if c==cLen
                                cPos = allPos;
                            else
                                cPos = find(cValsGood == cValsUnique(c));
                            end
                            
                            for t=1:tLen
                                if t==tLen
                                    tPos = allPos;
                                else
                                    tPos = find(tValsGood == tValsUnique(t));
                                end


                                aePos = intersect(aPos,ePos);
                                aesPos = intersect(aePos,sPos);
                                aesfPos = intersect(aesPos,fPos);
                                aesfoPos = intersect(aesfPos,oPos);
                                aesfocPos = intersect(aesfoPos,cPos);
                                aesfoctPos = intersect(aesfocPos,tPos);
                                parameterCombinations{a,e,s,f,o,c,t} = aesfoctPos; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end
    end

    % save
    save(fullfile(folderOut,'parameterCombinations.mat'),'parameters','parameterCombinations', ...
        'aValsUnique','eValsUnique','sValsUnique','fValsUnique','oValsUnique','cValsUnique','tValsUnique');

end
end
    
end

function [goodStimNums,goodStimTimes] = getGoodStimNumsGRF(folderOut,ignoreTargetStimFlag)

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