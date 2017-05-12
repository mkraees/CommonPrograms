% NEV stores information of only one of the Gabors (activeSide)
% 0 - Mapping 0
% 1 - Mapping 1
% 2 - Task Gabor

function matchingParameters=compareLLwithNEV(folderExtract,activeSide,showResults,convertToImageFlag)

load(fullfile(folderExtract,'LL.mat'));
load(fullfile(folderExtract,'stimResults.mat'));
load(fullfile(folderExtract,'digitalEvents.mat'));

% Compare basic properties
% NEV stores information of only one of the Gabors (activeSide)

if activeSide==2 % SRC protocol
    aziLL = LL.azimuthDeg;
    eleLL = LL.elevationDeg;
    sigmaLL = LL.sigmaDeg;
    radiusExists=0;
    sfLL = LL.spatialFreqCPD;
    oriLL = LL.orientationDeg;
    
elseif activeSide==0 % Map0
    if convertToImageFlag
        validMap = find(LL.stimType1==8);
    else
        validMap = find(LL.stimType1==1);
    end
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
    timeLL = LL.time1(validMap);
    
elseif activeSide==1 % Map2
    
    validMap = find(LL.stimType2==1);
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
    timeLL = LL.time2(validMap);
    
else % activeSide = [0 1]
    validMap = find(LL.stimType1==1);
    aziLL1 = LL.azimuthDeg1(validMap);
    eleLL1 = LL.elevationDeg1(validMap);
    sigmaLL1 = LL.sigmaDeg1(validMap);
    if isfield(LL,'radiusDeg1')
        radiusExists = 1;
        radiusLL1 = LL.radiusDeg1(validMap);
    else
        radiusExists = 0;
    end
    sfLL1 = LL.spatialFreqCPD1(validMap);
    oriLL1 = LL.orientationDeg1(validMap);
    conLL1 = LL.contrastPC1(validMap); 
    tfLL1 = LL.temporalFreqHz1(validMap);
    timeLL1 = LL.time1(validMap);
    
    validMap = find(LL.stimType2==1);
    aziLL2 = LL.azimuthDeg2(validMap);
    eleLL2 = LL.elevationDeg2(validMap);
    sigmaLL2 = LL.sigmaDeg2(validMap);
    if isfield(LL,'radiusDeg2')
        radiusExists = 1;
        radiusLL2 = LL.radiusDeg2(validMap);
    else
        radiusExists = 0;
    end
    sfLL2 = LL.spatialFreqCPD2(validMap);
    oriLL2 = LL.orientationDeg2(validMap);
    conLL2 = LL.contrastPC2(validMap); 
    tfLL2 = LL.temporalFreqHz2(validMap);
    timeLL2 = LL.time2(validMap);
    
    %----------
    % Combine all values
    % Combinations depend on the gabors which are not NULL, which in turn
    % depend on the specific protocol that is run.
    % For protocols where a gabor is hidden (i.e. null), this ignores its
    % parameter values and combines for the valid gabors only
    aziLL = [aziLL1; aziLL2]; aziLL = (reshape(aziLL,[],1))';
    eleLL = [eleLL1; eleLL2]; eleLL = (reshape(eleLL,[],1))';
    sigmaLL = [sigmaLL1; sigmaLL2]; sigmaLL = (reshape(sigmaLL,[],1))';
    sfLL = [sfLL1; sfLL2]; sfLL = (reshape(sfLL,[],1))';
    oriLL = [oriLL1; oriLL2]; oriLL = (reshape(oriLL,[],1))';
    conLL = [conLL1; conLL2]; conLL = (reshape(conLL,[],1))';
    tfLL = [tfLL1; tfLL2]; tfLL = (reshape(tfLL,[],1))';
    radiusLL = [radiusLL1; radiusLL2]; radiusLL = (reshape(radiusLL,[],1))';
    
    timeLL = timeLL1; % Vinay - if both gabors are drawn then assign timeLL1 (first drawn gabor) to timeLL
end

% Compare

if ~isfield(stimResults,'azimuth') %#ok<NODEF>
    disp('Stimulus parameter values not present in the digital data stream, taking from the LL file...');
    
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
    matchingParameters = [];
    
    % Saving updated stimResults 
    disp('Saving stimResults after taking values from LL file...');
    save(fullfile(folderExtract,'stimResults.mat'),'stimResults');
else

    if compareValues(aziLL,stimResults.azimuth)
        matchingParameters.azimuth=1;
        disp('Azimuths match.');
    else
        matchingParameters.azimuth=0; %#ok<*STRNU>
        error('*****************Azimuths do not match!!');
    end
    
    if compareValues(eleLL,stimResults.elevation)
        matchingParameters.elevation=1;
        disp('Elevations match.');
    else
        matchingParameters.elevation=0;
        disp('*****************Elevations do not match!!');
    end
    
    if compareValues(sigmaLL,stimResults.sigma)
        matchingParameters.sigma=1;
        disp('Sigmas match.');
    else
        matchingParameters.sigma=0;
        disp('*****************Sigmas do not match!!');
    end
    
    if radiusExists
        if compareValues(radiusLL,stimResults.radius)
            matchingParameters.radius=1;
            disp('Radii match.');
        else
            matchingParameters.radius=0;
            disp('******************Radii do not match!!');
        end
    else
        disp('Radius parameter does not exist');
    end
    
    if compareValues(sfLL,stimResults.spatialFrequency)
        matchingParameters.spatialFrequency=1;
        disp('Spatial Freq match.');
    else
        matchingParameters.spatialFrequency=0;
        disp('**************Spatial Freq do not match!!');
    end
    
    if compareValues(oriLL,stimResults.orientation)
        matchingParameters.orientation=1;
        disp('Orientations match.');
    else
        matchingParameters.orientation=0;
        disp('***************Orientations do not match!!');
    end
    
    if compareValues(conLL,stimResults.contrast)
        matchingParameters.contrast=1;
        disp('Contrasts match.');
    else
        matchingParameters.contrast=0;
        disp('*******************Contrasts do not match!!');
    end

    if compareValues(tfLL,stimResults.temporalFrequency)
        matchingParameters.temporalFrequency=1;
        disp('Temporal frequencies match.');
    else
        matchingParameters.temporalFrequency=0;
        disp('*******************Temporal frequencies do not match!!');
    end
end

if showResults
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Match start Times

    clear xD xL lxD lxL
    
    % Get NEV start times
    for i=1:length(digitalCodeInfo)
        if strcmp(digitalCodeInfo(i).codeName,'TS')
            endPos=i;
            break;
        end
    end

    if ~exist('endPos','var')
        error('No trialEvent named TS');
    else
        xD=digitalCodeInfo(endPos).time;
    end
    
    xL=LL.startTime;
    
    xD = diff(xD); lxD = length(xD); xD=xD(:);
    xL = diff(xL); lxL = length(xL); xL=xL(:);
    
    if lxD == lxL
        disp(['Number of startTrials: ' num2str(lxD)]);
        subplot(231)
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('Difference in Start Times (s)');
        legend('Dig','LL','Location','SouthEast');
        
        subplot(234)
        plot(1000*(xD-xL),'b');
        ylabel('Digital-Lablib times (ms)');
        xlabel('Trial Number');
        
        matchingParameters.maxChangeTrialsPercent = 100*max(abs(xD-xL) ./ xD);
    else
        disp(['Num of startTrials: digital: ' num2str(lxD+1) ' , LL: ' num2str(lxL+1)]);
        mlx = min(lxD,lxL);
        subplot(231)
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('Start Times (s)');
        legend('Dig','LL','Location','SouthEast');
        
        subplot(234)
        plot(1000*(xD(1:mlx)-xL(1:mlx)),'b');
        ylabel('Difference in start times (ms)');
        xlabel('Trial Number');
        
        matchingParameters.maxChangeTrialsPercent = 100*max(abs(xD(1:mlx)-xL(1:mlx)) ./ xD(1:mlx));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if activeSide~=2
        % Match stimulus Frame positions
        clear xD xL lxD lxL
        xD=stimResults.time;
        xL=timeLL/1000;
        
        xD = diff(xD); lxD = length(xD); xD=xD(:);
        xL = diff(xL); lxL = length(xL); xL=xL(:);
        
        if lxD == lxL
            disp(['Number of stimOnset: ' num2str(lxD)]);
            subplot(232)
            plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
            ylabel('Difference in StimOn time (s)');
            axis tight; 
            %legend('Dig','LL','Location','SouthEast');
            
            subplot(235)
            plot(1000*(xD-xL),'b');
            ylabel('Digital-Lablib times (ms)');
            xlabel('Stim Number');
            axis tight
            
            matchingParameters.maxChangeStimOnPercent = 100*max(abs(xD-xL) ./ xD);
        else
            disp(['Number of stimOnset: digital: ' num2str(lxD+1) ' , LL: ' num2str(lxL+1)]);
            mlx = min(lxD,lxL);
            subplot(232)
            plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
            ylabel('StartOn Frame');
            axis tight
            %legend('Dig','LL','Location','SouthEast');
            
            subplot(235)
            plot(1000*(xD(1:mlx)-xL(1:mlx)),'b');
            ylabel('Digital-Lablib times (ms)');
            xlabel('Stim Number');
            axis tight
            
            matchingParameters.maxChangeStimOnPercent = 100*max(abs(xD(1:mlx)-xL(1:mlx)) ./ xD(1:mlx));
        end
    end
    
    %%%% Match EOTCodes
    
    clear xD xL lxD lxL
    % Get NEV start times
    clear endPos
    for i=1:length(digitalCodeInfo)
        if strcmp(digitalCodeInfo(i).codeName,'TE')
            endPos=i;
            break;
        end
    end

    if ~exist('endPos','var')
        error('No trialEvent named TE');
    else
        xD=digitalCodeInfo(endPos).value;
    end
    xL=LL.eotCode;
    
    lxD = length(xD); xD=xD(:);
    lxL = length(xL); xL=xL(:);
    
    if lxD == lxL
        disp(['Number of eotCodes: ' num2str(lxD)]);
        subplot(233)
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('eotCode number');
        %legend('Dig','LL','Location','SouthEast');
        
        subplot(236)
        plot(xD-xL,'b');
        ylabel('\delta eotCode number');
        xlabel('Trial Number');
        
    else
        disp(['Number of stimOnset: digital: ' num2str(lxD) ' , LL: ' num2str(lxL)]);
        mlx = min(lxD,lxL);
        subplot(233)
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('eotCode number');
        axis tight
        %legend('Dig','LL','Location','SouthEast');
        
        subplot(236)
        plot(xD(1:mlx)-xL(1:mlx),'b');
        ylabel('\Delta eotCode number');
        xlabel('Trial Number');
        axis ([0 mlx -7 7]);
    end
end

end
function result = compareValues(x,y)
thres = 10^(-2);
if max(x-y) < thres;
    result=1;
else
    result=0;
end
end