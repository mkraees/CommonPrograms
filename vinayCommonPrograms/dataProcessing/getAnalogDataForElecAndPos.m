% get the analogData corresponding to any electrode and stimPos
% Vinay Shirhatti, 27 September 2016
%**************************************************************************

function elecData = getAnalogDataForElecAndPos(folderName,goodPos,electrodes,singleprecision,reftype,chanElecs,analogInputChannel)

if ~exist('reftype','var'); reftype = 'single'; end
if ~exist('singleprecision','var'); singleprecision = 0; end % By default use double precision. 
% Set this to 1 to use single precision
if ~exist('analogInputChannel','var'); analogInputChannel = 0; end

folderLFP = fullfile(folderName,'segmentedData','LFP');

elecData = cell(1,length(electrodes));
clear analogData analogDataX
for i = 1:length(electrodes)
    
    if iscell(goodPos)
        if length(goodPos)==1 % only common goodPos
            elecGoodPos = goodPos{1};
        else % electrodewise goodPos
            elecGoodPos = goodPos{i};
        end
    else
        elecGoodPos = goodPos;
    end
    
    if strcmpi(reftype,'single')
        if analogInputChannel
            load(fullfile(folderLFP,['ainp' num2str(electrodes(i)) '.mat']));
        else
            load(fullfile(folderLFP,['elec' num2str(electrodes(i)) '.mat']));
        end
        if ~singleprecision
            elecData{i} = analogData(elecGoodPos,:);
        else
            elecData{i} = single(analogData(elecGoodPos,:));% [Vinay] - seeing if single precision works, 27 Feb 2017
            % It works! :D
        end
    else
        load(fullfile(folderLFP,['elec' num2str(chanElecs(i,1)) '.mat']));
        analogDataX = analogData(elecGoodPos,:);
        clear analogData
        load(fullfile(folderLFP,['elec' num2str(chanElecs(i,2)) '.mat']));
        elecData{i} = analogDataX - analogData(elecGoodPos,:);
    end
end

end