% convert a cell to a row cell

function [out] = rowcell(in)

if size(in,1)~=1 && size(in,2)==1 && ndims(in)==2 %#ok<ISMAT>
    out = in';
else
    out = in;
end

end