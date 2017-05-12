% get all the valsUnique across the specified protocols
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function [a,e,s,f,o,c,t,r,p] = getCombinedValsUnique(folderSourceString,gridType,subjectslist,expdateslist,protocolslist)

a=[]; e=[]; s=[]; f=[]; o=[]; c=[]; t=[]; r=[]; p=[];

for i=1:length(subjectslist)
    folderName = fullfile(folderSourceString,'data',subjectslist{i},gridType,expdateslist{i},protocolslist{i});
    [~,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(fullfile(folderName,'extractedData'));
    a = cat(2,a,aValsUnique);
    e = cat(2,e,eValsUnique);
    s = cat(2,s,sValsUnique);
    f = cat(2,f,fValsUnique);
    o = cat(2,o,oValsUnique);
    c = cat(2,c,cValsUnique);
    t = cat(2,t,tValsUnique);
    r = cat(2,r,rValsUnique);
    p = cat(2,p,pValsUnique);
end

a = unique(a);
e = unique(e);
s = unique(s);
f = unique(f);
o = unique(o);
c = unique(c);
t = unique(t);
r = unique(r);
p = unique(p);

end