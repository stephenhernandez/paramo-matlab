function [templates,stnms,stla,stlo,ref,fs] = distillStruct(S,ncomps)
if nargin < 2
    ncomps = 3;
end

%%
ref = S(1).ref;
dt = S(1).delta;
fs = 1/dt;
if dt < 1
    fs = round(fs);
end

%%
lS = length(S);
ns = round(lS/ncomps);

%%
npts = S(1).npts;
stla = NaN(ns,1);
stlo = stla;
templates = NaN(npts,lS);
stnms = repmat("",ns,1);

%%
if ncomps == 1
    for i = 1:ns
        stla(i,1) = S(i).stla;
        stlo(i,1) = S(i).stlo;
        stnms(i,1) = S(i).kstnm;
        templates(:,i) = S(i).d;
        templates(:,i) = S(i).d;
        templates(:,i) = S(i).d;
    end
else
    for i = 1:ns
        stla(i,1) = S(1+(i-1)*3).stla;
        stlo(i,1) = S(1+(i-1)*3).stlo;
        stnms(i,1) = S(1+(i-1)*3).kstnm;
        templates(:,1+(i-1)*3) = S(1+(i-1)*3).d;
        templates(:,2+(i-1)*3) = S(2+(i-1)*3).d;
        templates(:,3+(i-1)*3) = S(3+(i-1)*3).d;
    end
end
