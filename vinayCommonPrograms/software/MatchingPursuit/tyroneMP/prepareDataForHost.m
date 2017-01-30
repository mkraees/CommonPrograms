% prepareDataForTyrone is a modified version of prepareDataForOrchestra 
% which is a generic program that converts the data in a
% format that can be run on the Tyrone cluster, SERC
%
% modified/cleaned from earlier prepareDataForHostVinay
% Vinay, 23 December 2015
%==========================================================================

function prepareDataForHost(monkeyName,expDate,protocolName,channelNumbers,folderSourceString,gridType,tRange,outputFolder)

folderNameMain = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName);

% Input folder
inputFolder  = fullfile(folderNameMain,'segmentedData','LFP');

if ~exist('tRange','var')
    tRange = [];
end

% Output folder
if ~exist('outputFolder','var')
    if ~isempty(tRange)
        outputFolder = fullfile(folderNameMain,'mpAnalysis','stimEpoch');
    else
        outputFolder = fullfile(folderNameMain,'mpAnalysis');
    end
    makeDirectory(outputFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Segmenting data...');
for i=1:length(channelNumbers)
    disp(channelNumbers(i));
    clear analogData analogInfo analogMatFile
    
    analogMatFile = fullfile(inputFolder,['elec' num2str(channelNumbers(i)) '.mat']);
    
    if ~exist(analogMatFile,'file')
        disp([analogMatFile ' does not exist']);
    else
        load(analogMatFile);
        X = analogData';
        if ~isempty(tRange)
            X = X(tRange,:);
        end

        if ~exist('Fs','var')
            Fs = analogInfo.SampleRate; % Initialize
        else
            if (Fs ~= analogInfo.SampleRate)
                error('Sampling rates not the same!!');
            end
        end
        
        if ~exist('L','var')
            L = size(X,1); % Initialize
        else
            if (size(X,1) ~= L)
                error('signal lengths not the same!!');
            end
        end

        if ispc
            tag = ['elec' num2str(channelNumbers(i)) '\'];
        else
            tag = ['elec' num2str(channelNumbers(i)) '/']; % Vinay - for linux
        end
        range = [1 L];
        importData(X,outputFolder,tag,range,Fs);

        Numb_points = L;
        Max_iterations = 500; % [Vinay] - changed from 500 to 100
        disp(['Preparing for electrode/channel ' num2str(channelNumbers(i))]);
        channelFolder = fullfile(outputFolder,['elec' num2str(channelNumbers(i))]);
        prepareMPForOrchestra(channelFolder,tag,Numb_points,Max_iterations,channelFolder);
    end
end

% write script file
% writeScriptFileTyrone(monkeyName,expDate,protocolName,channelNumbers,folderSourceString,gridType);
% script file not required for running this on the host machine
end