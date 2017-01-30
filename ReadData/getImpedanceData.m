function getImpedanceData(subjectName,expDate,folderSourceString,gridType)

folderNameRawData = fullfile(folderSourceString,'data','rawData',[subjectName expDate]);
folderNameSave = fullfile(folderSourceString,'data',subjectName,gridType,expDate);
makeDirectory(folderNameSave);

fileName = fullfile(folderNameRawData,[subjectName expDate 'Impedance']);
if ~exist(fileName,'file')    
    fileName = [fileName '.txt'];             
end

if ~exist(fileName,'file')
    disp([fileName ' does not exist']);
else
    X = textread(fileName,'%s'); %#ok<DTXTRD>
    
    impedanceValues = zeros(1,96);
    for i=1:96
        for j=1:length(X)
            if strcmp(X{j},['elec' num2str(i)])
                impedanceValues(i) = str2double(X{j+1});
                break;
            end
        end
    end
    % Save
    save(fullfile(folderNameSave,'impedanceValues.mat'),'impedanceValues');
end
end