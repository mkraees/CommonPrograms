% [Vinay] modified from showElectrodeLocationsInCustomColor
%
% ******* additional arguments:
% colorsList: list of colors to be used for the electrodes
% invertDefaultMap: flag to flip the blue and yellow corners of the default
% colormap
% colorFactor: lighten or darken all colors by this factor
%**************************************************************************
function [plotHandle,colorNames,electrodeArray,colorRows,colorCols] = showElectrodeLocationsInCustomColor(gridPosition,plotHandle,holdOnState,hideElectrodeNums,startRow,gridType,subjectName,gridLayout,colorsList,invertDefaultMap,colorFactor)

if ~exist('hideElectrodeNums','var');    hideElectrodeNums=0;           end
if ~exist('startRow','var');             startRow=5;                    end
if ~exist('gridType','var');             gridType = 'Microelectrode';   end
if ~exist('subjectName','var');          subjectName=[];                end
if ~exist('gridLayout','var');           gridLayout=2;                  end
if ~exist('colorsList','var');           
    colorsList=[];
    defColors = 1; % default colors
else
    if isempty(colorsList)
        defColors = 1;
    else
        defColors = 0;
    end
end
if ~exist('invertDefaultMap','var');     invertDefaultMap=0;            end
if ~exist('colorFactor','var');          colorFactor=1;                 end


if strcmpi(gridType,'ECoG')
    numRows=8;numCols=10;
elseif strcmpi(gridType,'Microelectrode')
    numRows=10;numCols=10;
else
    numRows=10;numCols=11;
end
[~,~,electrodeArray] = electrodePositionOnGrid(1,gridType,subjectName,gridLayout);

if ~exist('plotHandle','var') || isempty(plotHandle)
    plotHandle = subplot('Position',gridPosition,'XTickLabel',[],'YTickLabel',[],'box','on');
end
if ~exist('holdOnState','var')
    holdOnState = 1;
end

if ~holdOnState
    cla(plotHandle);
end

axes(plotHandle);
dX = 1/numCols;
dY = 1/numRows;

lineXRow = zeros(2,numRows);lineYRow = zeros(2,numRows);
for i=1:numRows
    lineXRow(:,i) = [0 1]; lineYRow(:,i) = [i*dY i*dY];
end
lineXCol = zeros(2,numCols);lineYCol = zeros(2,numCols);
for i=1:numCols
    lineXCol(:,i) = [i*dX i*dX]; lineYCol(:,i) = [0 1];
end
line(lineXRow,lineYRow,'color','k'); hold on
line(lineXCol,lineYCol,'color','k'); 
hold off;

% Color the patches - choose the subgrid to color
colorRows = startRow:numRows;
colorCols = 1:numCols;

% The four corners have the following colors

% yellow (1 1 0)                      % Green (0 1 0)

% red (1 0 0)                         % Blue (0 0 1)

numColorRows = length(colorRows);
numColorCols = length(colorCols);

colorNames = cell(numColorRows,numColorCols);
colorNum = 1;
for j=1:numColorRows
    for i=1:numColorCols
        
        row = colorRows(j);
        col = colorCols(i);
        
        electrodeNum = electrodeArray(row,col);
        
        if electrodeNum>0
            
            if defColors
                if invertDefaultMap
                    
                    % The four corners have the following colors

                    % blue (0 0 1)                      % Green (0 1 0)

                    % red (1 0 0)                         % yellow (1 1 0)
                    % red value does not depend on col, goes from 0 to 1
                    % with row
                    if numColorRows==1
                        rValue = 1-(i-1)/(numColorCols-1);
                    else
                        rValue = (j-1)/(numColorRows-1);
                    end

                    % green value goes from 0 to 1 with col, does not
                    % depend on row
                    if numColorCols==1
                        gValue = 0;
                    else
                        gValue = (i-1)/(numColorCols-1);
                    end

                    % blue value is 0 on the three corners and 1 at the top
                    % left corner. We approximate this by the following
                    if numColorRows==1
                        bNum1 = j; bDen1 = 1;
                    else
                        bNum1 = 1-(j-1)/numColorRows; bDen1 = 1;
                    end

                    if numColorCols==1
                        bNum2 = i; bDen2 = 1;
                    else
                        bNum2 = 1-(i-1)/numColorCols; bDen2 = 1;
                    end

                    bValue = sqrt((bNum1*bNum2)/(bDen1*bDen2));
                    colorName = [rValue gValue bValue];
                    if colorFactor
                        colorName = colorName*colorFactor;
                        colorName(colorName>1) = 1;
                    end

                else
                    % red value does not depend on row, goes from 1 to 0 with col
                    if numColorCols==1
                        rValue = 1;
                    else
                        rValue = 1-(i-1)/(numColorCols-1);
                    end

                    % green value goes from 1 to 0 with row, does not depend on col
                    if numColorRows==1
                        gValue = 0;
                    else
                        gValue = 1-(j-1)/(numColorRows-1);
                    end

                    % blue value is 0 on the three corners and 1 at the bottom right
                    % corner. We approximate this by the following
                    if numColorRows==1
                        bNum1 = j; bDen1 = 1;
                    else
                        bNum1 = j-1; bDen1 = numColorRows-1;
                    end

                    if numColorCols==1
                        bNum2 = i; bDen2 = 1;
                    else
                        bNum2 = i-1; bDen2 = numColorCols-1;
                    end

                    bValue = sqrt((bNum1*bNum2)/(bDen1*bDen2));
                    colorName = [rValue gValue bValue];
                    if colorFactor
                        colorName = colorName*colorFactor;
                        colorName(colorName>1) = 1;
                    end
                end
            else
                colorName = colorsList(colorNum,:);
                if colorFactor
                    colorName = colorName*colorFactor;
                    colorName(colorName>1) = 1;
                end
                colorNum = colorNum + 1;
            end
        else
            colorName = [1 1 1];
        end
            
        highlightRow=row;
        highlightCol=col;

        % Create patch
        patchX = (highlightCol-1)*dX;
        patchY = (numRows-highlightRow)*dY;
        patchLocX = [patchX patchX patchX+dX patchX+dX];
        patchLocY = [patchY patchY+dY patchY+dY patchY];
        
        patch(patchLocX,patchLocY,colorName);
        colorNames{row,col} = colorName;
    end
end

if ~hideElectrodeNums
    % Write electrode numbers
    for i=1:numRows
        textY = (numRows-i)*dY + dY/2;
        for j=1:numCols
            textX = (j-1)*dX + dX/2;
            if electrodeArray(i,j)>0
                text(textX,textY,num2str(electrodeArray(i,j)),'HorizontalAlignment','center');
            end
        end
    end
end
end