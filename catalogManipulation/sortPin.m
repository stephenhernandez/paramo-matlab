function [Psort,strdata] = sortPin(Pin)
%function [Psort,chardata,datadata] = sortPin(Pin)

%%
% if nargin < 2
%     sortFlag = true;
% end

%%
stnm = pull(Pin,'stnm');
ntwk = pull(Pin,'ntwk');
strdata = strcat(stnm,ntwk);
[~,sI] = sort(strdata);
Psort = Pin(sI);

%%
% Psort = Pin;
% nPphases = length(Pin);
% chardata = cell(nPphases,2);
% datadata = NaN(nPphases,8);
% 
% for i = 1:nPphases
%     chardata{i,1} = Pin(i).stnm;
%     chardata{i,2} = Pin(i).ntwk;
%     datadata(i,1) = Pin(i).dist;
%     datadata(i,2) = Pin(i).azimuth;
%     datadata(i,3) = Pin(i).hh;
%     datadata(i,4) = Pin(i).mmm;
%     datadata(i,5) = Pin(i).sec;
%     datadata(i,6) = Pin(i).res;
%     datadata(i,7) = Pin(i).wt;
%     datadata(i,8) = datenum(Pin(i).t);
% end
% 
% if sortFlag
%     [chardata,sI] = sortrows(chardata);
%     datadata = datadata(sI,:);
% end
% 
% for i = 1:nPphases
%     Psort(i).stnm = chardata{i,1};
%     Psort(i).ntwk = chardata{i,2};
%     Psort(i).dist = datadata(i,1);
%     Psort(i).azimuth = datadata(i,2);
%     Psort(i).hh = datadata(i,3);
%     Psort(i).mmm = datadata(i,4);
%     Psort(i).sec = datadata(i,5);
%     Psort(i).res = datadata(i,6);
%     Psort(i).wt = datadata(i,7);
%     Psort(i).t = dn2dt(datadata(i,8));
% end