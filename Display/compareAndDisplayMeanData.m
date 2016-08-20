% X - A cell array of length L, with each cell being a 2D array of size
% NixK (i=1...L). Ni is the number of stimuli. xs must be of length K.

% The format is similar to the program compareAndDisplayData. This prgoram,
% however, takes the average of the data between the time timeLims.

function p=compareAndDisplayMeanData(X,xs,hPlot,hSig,testMethod,timeLims,colorX,displayPlot)

defaultMethodForMoreThanTwoEntries   = 'anova';

if ~exist('hPlot','var');      hPlot = subplot(211);                    end
if ~exist('hSig','var');       hSig  = subplot(212);                    end
if ~exist('testMethod','var'); testMethod = 'default';                  end
if ~exist('displayPlot','var'); displayPlot=1;                          end

L=length(X);
if ~exist('colorX','var');     colorX = jet(L);                         end
if ~exist('xs','var');         xs = 1:size(X{1},2);                     end
if ~exist('timeLims','var');   timeLims = [xs(1) xs(end)];              end

% Pool data within the limits set by timeLims
goodPos = intersect(find(xs>=timeLims(1)),find(xs<timeLims(2)));

comparisonData=cell(1,L);
meanOfEachCondition=zeros(1,L);
for i=1:L
    clear tmpData
    tmpData = mean(X{i}(:,goodPos),2);
    comparisonData{i} = tmpData(:);
    meanOfEachCondition(i) = mean(tmpData);
end

%clear X

if L==2 
    if strcmpi(testMethod,'default') || strncmpi(testMethod,'ttest',5)
        disp('Significance computed using t-test');
        [~,p,~,stats] = ttest2(comparisonData{1},comparisonData{2},0.05,'both','unequal');
        
        testParameter = 'tstat';
        testParameterValue = num2str(stats.tstat);
        testDOF = num2str(stats.df);
        pValue  = num2str(p);
        
    else
        disp('Use testMethod ttest for faster and more accurate results');
        [p,table] = anova1ForCellArrays(comparisonData,testMethod);
        
        testParameter = table{1,5};
        testParameterValue = num2str(table{2,5});
        testDOF = [num2str(table{2,3})  ', ' num2str(table{3,3})];
        pValue  = num2str(p);
        
    end
else
    if strcmpi(testMethod,'default')
        [p,table] = anova1ForCellArrays(comparisonData,defaultMethodForMoreThanTwoEntries);
    elseif strncmpi(testMethod,'ttest',5)
        disp(['Can not use t-test for more than two entries. Using ' defaultMethodForMoreThanTwoEntries ' instead']);
        [p,table] = anova1ForCellArrays(comparisonData,defaultMethodForMoreThanTwoEntries);
    else
        [p,table] = anova1ForCellArrays(comparisonData,testMethod);
    end
    
    testParameter = table{1,5};
    testParameterValue = num2str(table{2,5});
    
    if strcmp(testParameter,'F')
        testDOF = [num2str(table{2,3})  ', ' num2str(table{3,3})];
    else
        testDOF = num2str(table{2,3});
    end
    pValue  = num2str(p);
    
end

if displayPlot
    % display data
    barPlotForCellArrays(comparisonData,1:L,colorX,hPlot,0);
    
    % display the significance results
    allData=[];
    for i=1:L
        allData=[allData ; comparisonData{i}]; %#ok<*AGROW>
        mCompariosonData(i) = mean(comparisonData{i});
    end
    
    tableEntries{1,1} = 'Mean';
    tableEntries{1,2} = num2str(mean(allData),3);
    disp(['means: ' num2str(mCompariosonData,3)]);
    
    tableEntries{2,1} = 'Max Diff';
    tableEntries{2,2} = num2str(max(meanOfEachCondition)-min(meanOfEachCondition),3);
    
    tableEntries{3,1} = 'Std';
    tableEntries{3,2} = num2str(std(allData),3);
    disp(['n for std = ' num2str(length(allData))]);
    
    tableEntries{4,1} = 'SEM';
    tableEntries{4,2} = num2str(std(allData)/sqrt(length(allData)));
    
    tableEntries{5,1} = testParameter;
    tableEntries{5,2} = testParameterValue;
    
    tableEntries{6,1} = 'DOF';
    tableEntries{6,2} = testDOF;
    
    tableEntries{7,1} = 'p';
    tableEntries{7,2} = pValue;
    
    makeTable(get(hSig,'Position'),tableEntries);
end
