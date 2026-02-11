function [P,t,meanIET,stdIET,Pamp] = pratio(t,nr,amps,meanFlag)
if nargin < 3
    amps = [];
end

if nargin < 4
    meanFlag = true;
end

%%
lT = length(t);

if ~isdatetime(t(1))
    t = dn2dt(t);
end
[IET,~,eI] = cutWindows(datenum(t),nr,nr-1,false);
IET = dn2dt(IET);
IET = seconds(diff(IET));
size(IET)

%%
%[IET,~,eI] = cutWindows(iet,nr,nr-1,false);

%%
if meanFlag
    meanIET = mean(IET);
    stdIET = std(IET,0);
else
    meanIET = median(IET);
    stdIET = mad(IET,1);
end

%%
t = t(eI);
P = (meanIET./stdIET)';

Pamp = [];
if ~isempty(amps)
    %check if size of amp and t vectors are the same
    lA = length(amps);
    if lT ~= lA
        return;
    end
    ampCut = cutWindows(amps,nr,nr-1,false);
    %ampCut = diff(ampCut);
    %size(ampCut)
    if meanFlag
        meanAmp = mean(ampCut);
        stdAmp = std(ampCut,0);
    else
        meanAmp = median(ampCut);
        stdAmp = mad(ampCut,1);
    end
    Pamp = (meanAmp./stdAmp)';
end
