function [Ssort,chardata,datadata] = sortSin(Sin,sortFlag)
if nargin < 2
    sortFlag = true;
end
Ssort = Sin;
nSphases = length(Sin);
chardata = cell(nSphases,2);
datadata = NaN(nSphases,8);

for i = 1:nSphases
    chardata{i,1} = Sin(i).stnm;
    chardata{i,2} = Sin(i).ntwk;
    datadata(i,1) = Sin(i).dist;
    datadata(i,2) = Sin(i).azimuth;
    datadata(i,3) = Sin(i).hh;
    datadata(i,4) = Sin(i).mmm;
    datadata(i,5) = Sin(i).sec;
    datadata(i,6) = Sin(i).res;
    datadata(i,7) = Sin(i).wt;
    datadata(i,8) = datenum(Sin(i).t);
end

if sortFlag
    [chardata,sI] = sortrows(chardata);
    datadata = datadata(sI,:);
end

for i = 1:nSphases
    Ssort(i).stnm = chardata{i,1};
    Ssort(i).ntwk = chardata{i,2};
    Ssort(i).dist = datadata(i,1);
    Ssort(i).azimuth = datadata(i,2);
    Ssort(i).hh = datadata(i,3);
    Ssort(i).mmm = datadata(i,4);
    Ssort(i).sec = datadata(i,5);
    Ssort(i).res = datadata(i,6);
    Ssort(i).wt = datadata(i,7);
    Ssort(i).t = dn2dt(datadata(i,8));
end