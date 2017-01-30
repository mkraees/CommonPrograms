% Generate Time-Frequency spectrum using Matching Pursuit
% Vinay Shirhatti, 17 March 2016
%
%**************************************************************************

function [mpSpectrum,elecMPSpectrum] = genMPSpectrum(folderName,channelNumbers,goodPos,L,numAtomsMP,mpSaveSpectrumName)
    

for i = 1:length(channelNumbers)
    
    clear gaborInfo
    
    mpFolder = fullfile(folderName,'mpAnalysis',['elec' num2str(channelNumbers(i))]);
    gaborInfoFile = fullfile(mpFolder,'gaborInfo.mat');
    if ~exist(gaborInfoFile,'file')
        error('Generate the gaborInfo files first !!! Use runMPDataGen');
    else
        load(gaborInfoFile);
    end
    
    if ~exist('goodPos','var')
        goodPos = 1:length(gaborInfo);
    end
    
    gaborInfoGoodPos = gaborInfo(goodPos);

    if ~exist('L','var')
        L = 4096;
    end
    wrap = [];
    atomList = (1:numAtomsMP);
    mpSpectrum = [];
    disp(['Reconstructing Energy from: ' num2str(numAtomsMP) ' atoms, and '  num2str(length(goodPos)) ' trials']);
    for m=1:length(goodPos)
        disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
        rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
        if m == 1
            mpSpectrum = rEnergy;
        else
            mpSpectrum = mpSpectrum + rEnergy;
        end
    end
    mpSpectrum = mpSpectrum/length(goodPos);
    
    elecMPSpectrum(:,:,i) = mpSpectrum;

    if exist('mpSaveSpectrumName','var')
        disp('saving spectrum data...');
        save(fullfile(mpFolder,mpSaveSpectrumName), 'mpSpectrum');
    end
    
    clear mpSpectrum
    
end

mpSpectrum = squeeze(mean(elecMPSpectrum,3));

if exist('mpSaveSpectrumName','var')
    disp('saving spectrum data...');
    meanSpectraFolder = fullfile(mpFolder,'meanSpectra');
    makeDirectory(meanSpectraFolder);
    save(fullfile(meanSpectraFolder,mpSaveSpectrumName), 'mpSpectrum');
end


end