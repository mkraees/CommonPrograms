% This program generates the parameterCombinations variable from the
% stimResults
function parameterCombinations = getDisplayCombinationsCRS(folderOut,goodStimNums)

load(fullfile(folderOut,'stimResults.mat'));

% Nine parameters are chosen: for CRS [Vinay]
% 1. Azimuth
% 2. Elevation
% 3. Sigma
% 4. Spatial Frequency
% 5. Orientation
% 6. Contrast
% 7. Temporal Frequency
% 8. Spatial Phase
% 9. Radius

% Parameters index
parameters{1} = 'azimuth';
parameters{2} = 'elevation';
parameters{3} = 'sigma';
parameters{4} = 'spatialFrequency';
parameters{5} = 'orientation';
parameters{6} = 'contrast';
parameters{7} = 'temporalFrequency'; %#ok<NASGU>
parameters{8} = 'spatialPhase'; % Vinay - for CRS
parameters{9} = 'radius'; % Vinay - for CRS
parameters{10} = 'kGaborNumber'; % Vinay - for CRS

% [Vinay] - read the parameter values for all the gabors together
aValsAllGabor  = stimResults.azimuth;
eValsAllGabor  = stimResults.elevation;
sValsAllGabor  = stimResults.sigma;
fValsAllGabor  = stimResults.spatialFrequency;
oValsAllGabor  = stimResults.orientation;
cValsAllGabor  = stimResults.contrast;
tValsAllGabor  = stimResults.temporalFrequency;
pValsAllGabor  = stimResults.spatialPhase; % Vinay - for CRS
rValsAllGabor  = stimResults.radius; % Vinay - for CRS

% [Vinay] - The parameters stored above are for all the gabors together.
% They appear in repeated ordered groups (gabor0,1,2). So the relevant
% parameter values have to be picked corresponding to individual gabors.

% [Vinay] - If there are any null gabors then assign '0' values to its
% parameters and all lengths as 1

numOfGabors = length(stimResults.side);

if numOfGabors ~= 3
    
    % Vinay - If the number of gabors is not 3 there are null gabors.
    % Reshape the ValsAllGabor matrices so that each row corresponds to the
    % values for a particular gabor.
    aTemp = reshape(aValsAllGabor,numOfGabors,length(aValsAllGabor)/numOfGabors);
    eTemp = reshape(eValsAllGabor,numOfGabors,length(eValsAllGabor)/numOfGabors);
    sTemp = reshape(sValsAllGabor,numOfGabors,length(sValsAllGabor)/numOfGabors);
    fTemp = reshape(fValsAllGabor,numOfGabors,length(fValsAllGabor)/numOfGabors);
    oTemp = reshape(oValsAllGabor,numOfGabors,length(oValsAllGabor)/numOfGabors);
    cTemp = reshape(cValsAllGabor,numOfGabors,length(cValsAllGabor)/numOfGabors);
    tTemp = reshape(tValsAllGabor,numOfGabors,length(tValsAllGabor)/numOfGabors);
    pTemp = reshape(pValsAllGabor,numOfGabors,length(pValsAllGabor)/numOfGabors);
    rTemp = reshape(rValsAllGabor,numOfGabors,length(rValsAllGabor)/numOfGabors);
    
    countGabor = 0;
    for gi = 1:3
        if isempty(stimResults.side(stimResults.side==(gi-1))) % check if this gabor is null or not
            % If it is null then assign '0's to the corresponding row in a
            % temp matrix below
            aValsAllGaborTmp(gi,1:size(aTemp,2))  = 0;
            eValsAllGaborTmp(gi,1:size(eTemp,2))  = 0;
            sValsAllGaborTmp(gi,1:size(sTemp,2))  = 0;
            fValsAllGaborTmp(gi,1:size(fTemp,2))  = 0;
            oValsAllGaborTmp(gi,1:size(oTemp,2))  = 0;
            cValsAllGaborTmp(gi,1:size(cTemp,2))  = 0;
            tValsAllGaborTmp(gi,1:size(tTemp,2))  = 0;
            pValsAllGaborTmp(gi,1:size(pTemp,2))  = 0;
            rValsAllGaborTmp(gi,1:size(rTemp,2))  = 0;
        else
            % If it is not null then extract the corresponding row from the
            % aTemp (or the related one) matrix and assign the values to
            % the corresponding row of the Temp matrix below
            countGabor=countGabor+1;
            aValsAllGaborTmp(gi,1:size(aTemp,2))  = aTemp(countGabor,1:size(aTemp,2));
            eValsAllGaborTmp(gi,1:size(eTemp,2))  = eTemp(countGabor,1:size(eTemp,2));
            sValsAllGaborTmp(gi,1:size(sTemp,2))  = sTemp(countGabor,1:size(sTemp,2));
            fValsAllGaborTmp(gi,1:size(fTemp,2))  = fTemp(countGabor,1:size(fTemp,2));
            oValsAllGaborTmp(gi,1:size(oTemp,2))  = oTemp(countGabor,1:size(oTemp,2));
            cValsAllGaborTmp(gi,1:size(cTemp,2))  = cTemp(countGabor,1:size(cTemp,2));
            tValsAllGaborTmp(gi,1:size(tTemp,2))  = tTemp(countGabor,1:size(tTemp,2));
            pValsAllGaborTmp(gi,1:size(pTemp,2))  = pTemp(countGabor,1:size(pTemp,2));
            rValsAllGaborTmp(gi,1:size(rTemp,2))  = rTemp(countGabor,1:size(rTemp,2));
        end
    end
    
    % Redefine the ValsAllGabor matrices by reshaping the Temp matrix 
    % obtained above 
    aValsAllGabor = reshape(aValsAllGaborTmp,1,[]);
    eValsAllGabor = reshape(eValsAllGaborTmp,1,[]);
    sValsAllGabor = reshape(sValsAllGaborTmp,1,[]);
    fValsAllGabor = reshape(fValsAllGaborTmp,1,[]);
    oValsAllGabor = reshape(oValsAllGaborTmp,1,[]);
    cValsAllGabor = reshape(cValsAllGaborTmp,1,[]);
    tValsAllGabor = reshape(tValsAllGaborTmp,1,[]);
    pValsAllGabor = reshape(pValsAllGaborTmp,1,[]);
    rValsAllGabor = reshape(rValsAllGaborTmp,1,[]);
    
end
    
for gaborNum = 1:3
    
    % [Vinay] - since the ValsAllGood matrices have been redefined to 
    % contain the information for all the gabors we increment by 3
    
    gaborIndices = gaborNum:3:length(aValsAllGabor); % [Vinay] - index values corresponding to 'gaborNum' gabor
    
    aValsAll = aValsAllGabor(gaborIndices);
    eValsAll = eValsAllGabor(gaborIndices);
    sValsAll = sValsAllGabor(gaborIndices);
    fValsAll = fValsAllGabor(gaborIndices);
    oValsAll = oValsAllGabor(gaborIndices);
    cValsAll = cValsAllGabor(gaborIndices);
    tValsAll = tValsAllGabor(gaborIndices);
    pValsAll = pValsAllGabor(gaborIndices); % Vinay - for CRS
    rValsAll = rValsAllGabor(gaborIndices); % Vinay - for CRS
    
    if ~isempty(aValsAll)
        % Get good stim
        if ~exist('goodStimNums','var')
            goodStimNums = getGoodStimNumsCRS(folderOut);
        end

        aValsGood = aValsAll(goodStimNums);
        eValsGood = eValsAll(goodStimNums);
        sValsGood = sValsAll(goodStimNums);
        fValsGood = fValsAll(goodStimNums);
        oValsGood = oValsAll(goodStimNums);
        cValsGood = cValsAll(goodStimNums);
        tValsGood = tValsAll(goodStimNums);
        pValsGood = pValsAll(goodStimNums); % Vinay - for CRS
        rValsGood = rValsAll(goodStimNums); % Vinay - for CRS

        aValsUnique{gaborNum} = unique(aValsGood); aLen = length(aValsUnique{gaborNum});
        eValsUnique{gaborNum} = unique(eValsGood); eLen = length(eValsUnique{gaborNum});
        sValsUnique{gaborNum} = unique(sValsGood); sLen = length(sValsUnique{gaborNum});
        fValsUnique{gaborNum} = unique(fValsGood); fLen = length(fValsUnique{gaborNum});
        oValsUnique{gaborNum} = unique(oValsGood); oLen = length(oValsUnique{gaborNum});
        cValsUnique{gaborNum} = unique(cValsGood); cLen = length(cValsUnique{gaborNum});
        tValsUnique{gaborNum} = unique(tValsGood); tLen = length(tValsUnique{gaborNum});
        pValsUnique{gaborNum} = unique(pValsGood); pLen = length(pValsUnique{gaborNum}); % Vinay - for CRS
        rValsUnique{gaborNum} = unique(rValsGood); rLen = length(rValsUnique{gaborNum}); % Vinay - for CRS

        % display
        disp(['Gabor ' num2str(gaborNum - 1)]); % Vinay - display kGabor number
        disp(['Number of unique azimuths: ' num2str(aLen)]);
        disp(['Number of unique elevations: ' num2str(eLen)]);
        disp(['Number of unique sigmas: ' num2str(sLen)]);
        disp(['Number of unique Spatial freqs: ' num2str(fLen)]);
        disp(['Number of unique orientations: ' num2str(oLen)]);
        disp(['Number of unique contrasts: ' num2str(cLen)]);
        disp(['Number of unique temporal freqs: ' num2str(tLen)]);
        disp(['Number of unique spatial phases: ' num2str(pLen)]); % Vinay - for CRS
        disp(['Number of unique radii: ' num2str(rLen)]); % Vinay - for CRS

        % If more than one value, make another entry with all values
        if (aLen > 1)           aLen=aLen+1;                    end
        if (eLen > 1)           eLen=eLen+1;                    end
        if (sLen > 1)           sLen=sLen+1;                    end
        if (fLen > 1)           fLen=fLen+1;                    end
        if (oLen > 1)           oLen=oLen+1;                    end
        if (cLen > 1)           cLen=cLen+1;                    end
        if (tLen > 1)           tLen=tLen+1;                    end
        if (pLen > 1)           pLen=pLen+1;                    end % Vinay - for CRS
        if (rLen > 1)           rLen=rLen+1;                    end % Vinay - for CRS

        allPos = 1:length(goodStimNums);
        disp(['total combinations: ' num2str((aLen)*(eLen)*(sLen)*(fLen)*(oLen)*(cLen)*(tLen)*(pLen)*(rLen))]);

        for a=1:aLen
            if a==aLen
                aPos = allPos;
            else
                aPos = find(aValsGood == aValsUnique{gaborNum}(a));
            end

            for e=1:eLen
                if e==eLen
                    ePos = allPos;
                else
                    ePos = find(eValsGood == eValsUnique{gaborNum}(e));
                end

                for s=1:sLen
                    if s==sLen
                        sPos = allPos;
                    else
                        sPos = find(sValsGood == sValsUnique{gaborNum}(s));
                    end

                    for f=1:fLen
                        if f==fLen
                            fPos = allPos;
                        else
                            fPos = find(fValsGood == fValsUnique{gaborNum}(f));
                        end

                        for o=1:oLen
                            if o==oLen
                                oPos = allPos;
                            else
                                oPos = find(oValsGood == oValsUnique{gaborNum}(o));
                            end

                            for c=1:cLen
                                if c==cLen
                                    cPos = allPos;
                                else
                                    cPos = find(cValsGood == cValsUnique{gaborNum}(c));
                                end

                                for t=1:tLen
                                    if t==tLen
                                        tPos = allPos;
                                    else
                                        tPos = find(tValsGood == tValsUnique{gaborNum}(t));
                                    end

                                    for p=1:pLen % Vinay - for CRS
                                        if p==pLen
                                            pPos = allPos;
                                        else
                                            pPos = find(pValsGood == pValsUnique{gaborNum}(p));
                                        end

                                        for r=1:rLen % Vinay - for CRS
                                            if r==rLen
                                                rPos = allPos;
                                            else
                                                rPos = find(rValsGood == rValsUnique{gaborNum}(r));
                                            end


                                    aePos = intersect(aPos,ePos);
                                    aesPos = intersect(aePos,sPos);
                                    aesfPos = intersect(aesPos,fPos);
                                    aesfoPos = intersect(aesfPos,oPos);
                                    aesfocPos = intersect(aesfoPos,cPos);
                                    aesfoctPos = intersect(aesfocPos,tPos);
                                    aesfoctpPos = intersect(aesfoctPos,pPos); % Vinay - for CRS
                                    aesfoctprPos = intersect(aesfoctpPos,rPos); % Vinay - for CRS

                                    parameterCombinations{a,e,s,f,o,c,t,p,r,gaborNum} = aesfoctprPos; %#ok<AGROW> % [Vinay] - added additional parameters p and r for CRS
                                    % [Vinay] - the last dimension has been
                                    % added to save the gabor number
                                    % corresponding to its parameter
                                    % combinations. gabor0 (surround) has index
                                    % value 1, gabor1 (ring) -> 2,
                                    % gabor2(centre) -> 3
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

end

    % save
    save(fullfile(folderOut,'parameterCombinations.mat'),'parameters','parameterCombinations', ...
        'aValsUnique','eValsUnique','sValsUnique','fValsUnique','oValsUnique','cValsUnique','tValsUnique','pValsUnique','rValsUnique'); % [Vinay] - added additional parameters p and r for CRS

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