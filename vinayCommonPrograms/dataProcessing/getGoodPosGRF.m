% get goodPos for any protocol
% Vinay Shirhatti, 26 September 2016
% modified from a local function inside figure1HueTuningv3
% *************************************************************************

function [goodPos, badTrials] = getGoodPosGRF(folderName,vals,removeCommonBadtrials,imageStim,imagesParamList)

if ~exist('imagesParamList','var')
    imagesParamList = [];
end

if ~exist('imageStim','var')
    imageStim = 0;
end

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);

a = vals{1};
e = vals{2};
s = vals{3};
f = vals{4};
o = vals{5};
c = vals{6};
t = vals{7};

%**************************************************************************
% for every parameter, check if a value is passed. If it is empty then
% just take the collection of all values i.e. all pValsUnique
if isempty(a) 
    aPos = 1:length(aValsUnique);
else
    aPos = find(ismember(aValsUnique,a));
end

if isempty(e)
    ePos=1:length(eValsUnique);
else
    ePos = find(ismember(eValsUnique,e));
end

if isempty(s)
    sPos = 1:length(sValsUnique);
else
    if sum(ismember(sValsUnique,s))
        sPos = find(ismember(sValsUnique,s));
    else
        if sum(ismember(round(sValsUnique,1),round(s,1)))
            sPos = find(ismember(round(sValsUnique,1),round(s,1)));
        end
    end
end

if isempty(t)
    tPos = length(tValsUnique);
else
    tPos = find(ismember(tValsUnique,t));
end

if imageStim
    if isempty(o)
        oPos = 1:length(oValsUnique);
    else
        oPos = find(ismember(imagesParamList(1,:),o));
    end
    
    if isempty(f)
        fPos = 1:length(fValsUnique);
    else
        fPos = find(ismember(imagesParamList(2,:),f));
    end
    
    if isempty(c)
        cPos = 1:length(cValsUnique);
    else
        cPos = find(ismember(imagesParamList(3,:),c));
    end
    
    fPos = intersect(oPos,fPos); % imageNumber is mapped to sf
    fPos = intersect(cPos,fPos); % imageNumber is mapped to sf
else
    if isempty(f)
        fPos = 1:length(fValsUnique);
    else
        fPos = find(ismember(fValsUnique,f));
    end
    
    if isempty(o)
        oPos = 1:length(oValsUnique);
    else
        oPos = find(ismember(oValsUnique,o));
    end
    
    if isempty(c)
        cPos = 1:length(cValsUnique);
    else
        cPos = find(ismember(cValsUnique,c));
    end
end

clear goodPos

goodPos=[];

if ~isempty(aPos) && ~isempty(ePos) && ~isempty(sPos) && ~isempty(fPos) && ~isempty(oPos) && ~isempty(cPos) && ~isempty(tPos)
    goodPos = [parameterCombinations{aPos,ePos,sPos,fPos,oPos,cPos,tPos}]; 
    % the [] brackets combine all the goodPos from multiple Pos for any 
    % parameter here (which happens in cases of pooling across values)
end

if removeCommonBadtrials
    % Get bad trials
    badTrialFile = fullfile(folderSegment,'badTrials.mat');
    if ~exist(badTrialFile,'file')
        disp('Bad trial file does not exist...');
        badTrials=[];
    else
        [~, badTrials] = loadBadTrials(badTrialFile);
        disp([num2str(length(badTrials)) ' bad trials']);
    end

    goodPos = setdiff(goodPos,badTrials);
end
end
%==========================================================================
function [allBadTrials, badTrials, nameElec] = loadBadTrials(badTrialFile)
load(badTrialFile);
end