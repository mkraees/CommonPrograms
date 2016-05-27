function plotHandle = showElectrodeLocations(gridPosition,highlightElectrodes,colorNames,plotHandle,holdOnState,hideElectrodeNums,gridType,subjectName,gridLayout)

if ~exist('subjectName','var');          subjectName=[];                end
if ~exist('hideElectrodeNums','var');    hideElectrodeNums=0;           end
if ~exist('gridType','var');             gridType = 'Microelectrode';   end
if ~exist('gridLayout','var');           gridLayout=2;                  end

[~,~,electrodeArray] = electrodePositionOnGrid(1,gridType,subjectName,gridLayout);
[numRows,numCols] = size(electrodeArray);

if ~exist('plotHandle','var') || isempty(plotHandle)
    plotHandle = subplot('Position',gridPosition,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[],'box','on');
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
line(lineXRow,lineYRow,'color','k'); hold on;
line(lineXCol,lineYCol,'color','k'); 
hold off;

if ~isempty(highlightElectrodes)
    for i=1:length(highlightElectrodes)
        highlightElectrode=highlightElectrodes(i);
        
        [highlightRow,highlightCol,electrodeArray] = electrodePositionOnGrid(highlightElectrode,gridType,subjectName,gridLayout);

        % Create patch
        patchX = (highlightCol-1)*dX;
        patchY = (numRows-highlightRow)*dY;
        patchLocX = [patchX patchX patchX+dX patchX+dX];
        patchLocY = [patchY patchY+dY patchY+dY patchY];
        
        if iscell(colorNames)
            patch('XData',patchLocX,'YData',patchLocY,'FaceColor',colorNames{i});
        else
            patch('XData',patchLocX,'YData',patchLocY,'FaceColor',colorNames);
        end
    end
else
    [~,~,electrodeArray] = electrodePositionOnGrid(1,gridType,subjectName,gridLayout);
end

% Write electrode numbers
if ~hideElectrodeNums
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
set(plotHandle,'XTickLabel',[],'YTickLabel',[]);