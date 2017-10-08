% draw color stimuli as per the specified parameters
% taken out from computAndPlotMeanERPMeasures and made into a separate
% function, Vinay Shirhatti, 22 October 2016
%
% stimulus parameters are read from params.valsMatrix (which stores the
% values for every stimulus parameter for every condition) and the stimuli
% are drawn for every condition accordingly
%==========================================================================

function [params,stimColorPatch] = drawColorStimuliFromValsMatrix(params,numRows,numCols)

if ~exist('numRows','var')
    numRows = size(params.valsMatrix,1);
end

if ~exist('numCols','var')
    numCols = size(params.valsMatrix,2);
end

% check if the plot handles exist, otherwise create a new figure with the
% plot handles
if ~isfield(params,'hStimFig')
    figure;
    plotsPos = [0.08 0.15 0.87 0.8];
    if numCols>params.maxPlotsAlongX
        cols = params.maxPlotsAlongX; rows = numRows*ceil(numCols/cols);
        subplotsRearrangement = 1;
    else
        cols = numCols; rows = numRows; subplotsRearrangement = 0;
    end
    params.hStimFig = getPlotHandles(rows,cols,plotsPos);
else
    cols = numCols; rows = numRows; subplotsRearrangement = 0;
end

stimColorPatch = cell(numRows,numCols);
for r = 1:numRows
    for c = 1:numCols
         if ~subplotsRearrangement
            x = r; y = c;
            subplot(params.hStimFig(x,y));
        else
            x = (r-1)*ceil(numCols/cols)+ceil(c/cols);
            y = mod(c,cols);
            if y==0; y = cols; end
            subplot(params.hStimFig(x,y));
        end
        stim.azimuthDeg = params.valsMatrix{r,c,1}; % a
        stim.elevationDeg = params.valsMatrix{r,c,2}; % e
        stim.spatialFreqCPD = params.valsMatrix{r,c,4}; % saturation
        stim.orientationDeg = params.valsMatrix{r,c,5}; % hue
        stim.contrastPC = params.valsMatrix{r,c,6}; % value
        stim.radiusDeg = params.valsMatrix{r,c,3}*3; % radius = sigma*3

        stimColorPatch{r,c} = makeColourPatchHSV(params.hStimFig(x,y),stim);

    end
end

% Set the xlims and ylims as per the monitor diensions and viewing
% distance
%----------------------------------------------------------------------
% Dimensions and distances
%----------------------------------------------------------------------

distSub = 21.3; % subject's distance from the screen: 54 cm

% Monitor dimensions in inches
ht = 11.8; 
wd = 20.9;

% Resolution (num of pixels)
xRes = 1280;
yRes = 720;

%----------------------------------------------------------------------
% Calculations
%----------------------------------------------------------------------

yDeg = 2*(atan((ht/2)/distSub))*180/pi; % in deg
xDeg = 2*(atan((wd/2)/distSub))*180/pi;

% set the xlims and ylims as per monitor dimensions
set(params.hStimFig,'xlim',[-xDeg/2 xDeg/2]);
set(params.hStimFig,'ylim',[-yDeg/2 yDeg/2]);

% delete unused axes
for r=1:rows
    for c = 1:cols
        if cols*(r-1)+c>numRows*numCols
            axisHandle = findall(params.hStimFig(r,c));
            delete(axisHandle);
        end
    end
end

end
