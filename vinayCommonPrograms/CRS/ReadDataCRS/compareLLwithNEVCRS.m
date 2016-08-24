% NEV stores information of only one of the Gabors (activeSide)
% 0 - Mapping 0
% 1 - Mapping 1
% 2 - Task Gabor
% *** The above is for GRF. For CRS:
% 0 - Task gabor
% 1 - gabor0 (S) (surround gabor)
% 2 - gabor1 (R) (ring gabor)
% 3 - gabor2 (C) (centre gabor)

function matchingParameters=compareLLwithNEVCRS(folderExtract,~,showResults)

load(fullfile(folderExtract,'LL.mat'));
load(fullfile(folderExtract,'StimResults.mat'));
load(fullfile(folderExtract,'digitalEvents.mat'));

% Compare basic properties

% [Vinay] - In CRS index 0 in LL represents the task gabor, whereas index
% 1, 2 and 3 represent gabor0 (S), gabor1 (R) and gabor2 (C) respectively
% activeSide can be 0,1,2,[0 1],[0 2],[1 2] or [0 1 2]
% Irrespective of that, the parameter values for each gabor are stored.
% Therefore we just get the values for each one of them.
    
    validMap = find(LL.stimType1==1);
    aziLL1 = LL.azimuthDeg1(validMap);
    eleLL1 = LL.elevationDeg1(validMap);
    sigmaLL1 = LL.sigmaDeg1(validMap);
    radiusLL1 = LL.radiusDeg1(validMap);
    sfLL1 = LL.spatialFreqCPD1(validMap);
    oriLL1 = LL.orientationDeg1(validMap);
    conLL1 = LL.contrastPC1(validMap); 
    tfLL1 = LL.temporalFreqHz1(validMap);
    spLL1 = LL.SpatialPhaseDeg1(validMap);
    timeLL1 = LL.time1(validMap);
    
    validMap = find(LL.stimType2==1);
    aziLL2 = LL.azimuthDeg2(validMap);
    eleLL2 = LL.elevationDeg2(validMap);
    sigmaLL2 = LL.sigmaDeg2(validMap);
    radiusLL2 = LL.radiusDeg2(validMap);
    sfLL2 = LL.spatialFreqCPD2(validMap);
    oriLL2 = LL.orientationDeg2(validMap);
    conLL2 = LL.contrastPC2(validMap); 
    tfLL2 = LL.temporalFreqHz2(validMap);
    spLL2 = LL.SpatialPhaseDeg2(validMap);
    timeLL2 = LL.time2(validMap);
    
    validMap = find(LL.stimType3==1);
    aziLL3 = LL.azimuthDeg3(validMap);
    eleLL3 = LL.elevationDeg3(validMap);
    sigmaLL3 = LL.sigmaDeg3(validMap);
    radiusLL3 = LL.radiusDeg3(validMap);
    sfLL3 = LL.spatialFreqCPD3(validMap);
    oriLL3 = LL.orientationDeg3(validMap);
    conLL3 = LL.contrastPC3(validMap); 
    tfLL3 = LL.temporalFreqHz3(validMap);
    spLL3 = LL.SpatialPhaseDeg3(validMap);
    timeLL3 = LL.time3(validMap);
    
    %----------
    % Combine all values
    % Combinations depend on the gabors which are not NULL, which in turn
    % depend on the specific protocol that is run.
    % For protocols where a gabor is hidden (i.e. null), this ignores its
    % parameter values and combines for the valid gabors only
    aziLL = [aziLL1; aziLL2; aziLL3]; aziLL = (reshape(aziLL,[],1))';
    eleLL = [eleLL1; eleLL2; eleLL3]; eleLL = (reshape(eleLL,[],1))';
    sigmaLL = [sigmaLL1; sigmaLL2; sigmaLL3]; sigmaLL = (reshape(sigmaLL,[],1))';
    sfLL = [sfLL1; sfLL2; sfLL3]; sfLL = (reshape(sfLL,[],1))';
    oriLL = [oriLL1; oriLL2; oriLL3]; oriLL = (reshape(oriLL,[],1))';
    conLL = [conLL1; conLL2; conLL3]; conLL = (reshape(conLL,[],1))';
    tfLL = [tfLL1; tfLL2; tfLL3]; tfLL = (reshape(tfLL,[],1))';
    spLL = [spLL1; spLL2; spLL3]; spLL = (reshape(spLL,[],1))';
    radiusLL = [radiusLL1; radiusLL2; radiusLL3]; radiusLL = (reshape(radiusLL,[],1))';
    
% Compare

if ~isfield(stimResults,'azimuth') %#ok<NODEF>
    disp('Stimulus parameter values not present in the digital data stream, taking from the LL file...');
    
    stimResults.azimuth = aziLL;
    stimResults.elevation = eleLL;
    stimResults.contrast = conLL;
    stimResults.temporalFrequency = tfLL;
    stimResults.radius = radiusLL;
    stimResults.sigma = sigmaLL;
    stimResults.orientation = oriLL;
    stimResults.spatialFrequency = sfLL;
    stimResults.spatialPhase = spLL;
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
    
    if compareValues(radiusLL,stimResults.radius)
        matchingParameters.radius=1;
        disp('Radii match.');
    else
        matchingParameters.radius=0;
        disp('******************Radii do not match!!');
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
    
    if compareValues(spLL,stimResults.spatialPhase)
        matchingParameters.spatialPhase=1;
        disp('Spatial phases match.');
    else
        matchingParameters.spatialPhase=0;
        disp('***************Spatial phases do not match!!');
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
        subplot(4,6,[1 2 7 8])
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('Difference in Start Times (s)');
        legend('Dig','LL','Location','SouthEast');
        
        subplot(4,6,[13 14 19 20])
        plot(1000*(xD-xL),'b');
        ylabel('Digital-Lablib times (ms)');
        xlabel('Trial Number');
        
        matchingParameters.maxChangeTrialsPercent = 100*max(abs(xD-xL) ./ xD);
    else
        disp(['Num of startTrials: digital: ' num2str(lxD+1) ' , LL: ' num2str(lxL+1)]);
        mlx = min(lxD,lxL);
        subplot(4,6,[1 2 7 8])
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('Start Times (s)');
        legend('Dig','LL','Location','SouthEast');
        
        subplot(4,6,[13 14 19 20])
        plot(1000*(xD(1:mlx)-xL(1:mlx)),'b');
        ylabel('Difference in start times (ms)');
        xlabel('Trial Number');
        
        matchingParameters.maxChangeTrialsPercent = 100*max(abs(xD(1:mlx)-xL(1:mlx)) ./ xD(1:mlx));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for gi = 1:3
        % Match stimulus Frame positions
        clear xD xL lxD lxL
        
        if gi == 1
            xD=stimResults.time0;
            xL=timeLL1/1000;
        elseif gi == 2
            xD=stimResults.time1;
            xL=timeLL2/1000;
        elseif gi == 3
            xD=stimResults.time2;
            xL=timeLL3/1000;
        end
        
        xD = diff(xD); lxD = length(xD); xD=xD(:);
        xL = diff(xL); lxL = length(xL); xL=xL(:);
        
        if lxD == lxL
            disp(['Number of stimOnset: ' num2str(lxD)]);
            if gi == 1
                subplot(4,6,3);
            elseif gi == 2
                subplot(4,6,4);
            elseif gi == 3
                subplot(4,6,9);
            end
            plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
            ylabel('Difference in StimOn time (s)');
            axis tight; 
            %legend('Dig','LL','Location','SouthEast');
            
            if gi == 1
                subplot(4,6,15);
            elseif gi == 2
                subplot(4,6,16);
            elseif gi == 3
                subplot(4,6,21);
            end
            plot(1000*(xD-xL),'b');
            ylabel('Digital-Lablib times (ms)');
            xlabel('Stim Number');
            axis tight
            
            matchingParameters.maxChangeStimOnPercent = 100*max(abs(xD-xL) ./ xD);
        else
            disp(['Number of stimOnset: digital: ' num2str(lxD+1) ' , LL: ' num2str(lxL+1)]);
            mlx = min(lxD,lxL);
            if gi == 1
                subplot(4,6,3);
            elseif gi == 2
                subplot(4,6,4);
            elseif gi == 3
                subplot(4,6,9);
            end
            plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
            ylabel('StartOn Frame');
            axis tight
            %legend('Dig','LL','Location','SouthEast');
            
            if gi == 1
                subplot(4,6,15);
            elseif gi == 2
                subplot(4,6,16);
            elseif gi == 3
                subplot(4,6,21);
            end
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
        subplot(4,6,[5 6 11 12])
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('eotCode number');
        %legend('Dig','LL','Location','SouthEast');
        
        subplot(4,6,[17 18 23 24])
        plot(xD-xL,'b');
        ylabel('\delta eotCode number');
        xlabel('Trial Number');
        
    else
        disp(['Number of stimOnset: digital: ' num2str(lxD) ' , LL: ' num2str(lxL)]);
        mlx = min(lxD,lxL);
        subplot(4,6,[5 6 11 12])
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('eotCode number');
        axis tight
        %legend('Dig','LL','Location','SouthEast');
        
        subplot(4,6,[17 18 23 24])
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