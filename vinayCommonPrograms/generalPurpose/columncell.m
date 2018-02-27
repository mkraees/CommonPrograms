% convert a cell to a column cell

function [out] = columncell(in)

if size(in,2)~=1 && size(in,1)==1
    out = in';
else
    out = in;
end

end