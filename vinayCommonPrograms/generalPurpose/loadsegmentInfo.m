% load spikes segments info for a protocol
% modified from loadspikesInfo.m
% Vinay Shirhatti, 17 Jan 2017
%**************************************************************************

function [segmentChannelsStored,numItems] = loadsegmentInfo(folderSegments)
fileName = fullfile(folderSegments,'segmentInfo.mat');
if exist(fileName,'file')
    load(fileName);
else
    segmentChannelsStored=[];
    numItems=0;
end
end