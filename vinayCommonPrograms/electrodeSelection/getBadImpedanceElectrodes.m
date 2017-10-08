% get the list of electrodes with bad impedances
% taken from Poojya's figure codes for DualGammaPaper
% December 2016
% 
% added option to pass impedance cut-off limits, 29 April 2017
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function badElecs = getBadImpedanceElectrodes(folderSourceString,subjectName,expDate,gridType,impedanceCutoffs)

if ~exist('gridType','var'); gridType = 'Microelectrode'; end
if ~exist('impedanceCutoffs','var'); impedanceCutoffs = [250 2500]; end

try
    load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,'impedanceValues.mat'),'impedanceValues');
    badElecsList1=find(impedanceValues>impedanceCutoffs(2));
    badElecsList2=find(impedanceValues<impedanceCutoffs(1));
    badElecs=[badElecsList1 badElecsList2];
catch 
    if strcmpi(subjectName,'alpa')
       disp(' no impedance data available for this monkey');
       badElecs=[];
    else %% for kesari impedance data is available for all the days in which the datasets were collected
        disp('Saving Impedance data first'); 
        % make sure the commonPrograms used in the lab is mapped into the  path
        try
            getImpedanceData(subjectName,expDate,folderSourceString,gridType);
            load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,'impedanceValues.mat'),'impedanceValues');
            badElecsList1=find(impedanceValues>impedanceCutoffs(2));
            badElecsList2=find(impedanceValues<impedanceCutoffs(1));
            badElecs=[badElecsList1 badElecsList2];
        catch
            disp(' no impedance data available for this day');
            badElecs=[];
        end
    end
end 

end
