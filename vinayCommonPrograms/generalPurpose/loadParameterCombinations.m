% load parameter combinations for a protocol
% This function is the same as the one in display programs
% Vinay Shirhatti, 01 February 2016
%**************************************************************************


function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique, rValsUnique, pValsUnique] = loadParameterCombinations(folderExtract)
% [Vinay] - added rValsUnique and pValsUnique for radius and spatial phase
load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('rValsUnique','var')
    rValsUnique=[];
end

if ~exist('pValsUnique','var')
    pValsUnique=[];
end

if ~exist('sValsUnique','var')
    sValsUnique=rValsUnique/3;
end

if ~exist('cValsUnique','var')
    cValsUnique=[];
end

if ~exist('tValsUnique','var')
    tValsUnique=[];
end
end

