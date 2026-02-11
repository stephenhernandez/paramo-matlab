function S = shiftWaveforms(S,pTime,threeFlag)
%
% shiftWaveforms shift waveforms accoring to vector of arrival times
%
% S = shiftWaveforms(S,pTime,threeFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
if nargin < 3
    threeFlag = false;
end

%%
if ~isdatetime(pTime)
    disp('arrival times should be datetime structures');
    return;
end

%sizeS = size(S);
S = S(:);
lS = length(S);
lp = length(pTime);
minPtime = min(pTime);

%%
if threeFlag
    if lS ~= 3*lp
        disp('dimension mismatch');
        return;
    end
elseif lS ~= lp
    disp('dimension mismatch');
    return;
end

%%
delta = pull(S,'delta');
if any(diff(delta))
    if threeFlag
        disp('not all sample rates are the same');
        return
    end
    shift = zeros(lS,1);
    for i = 1:lS
        shift(i) = t2i(pTime(i),minPtime,delta(i));
    end
else
    delta = delta(1);
    shift = t2i(pTime,minPtime,delta);
end

%%
if threeFlag
    lS = lS/3;
    for i = 1:lS
        dtmp = S(3*(i-1)+1).d;
        dtmp = circshift(dtmp,-(shift(i)-1));
        S(3*(i-1)+1).d = dtmp;
        
        dtmp = S(3*(i-1)+2).d;
        dtmp = circshift(dtmp,-(shift(i)-1));
        S(3*(i-1)+2).d = dtmp;
        
        dtmp = S(3*(i-1)+3).d;
        dtmp = circshift(dtmp,-(shift(i)-1));
        S(3*(i-1)+3).d = dtmp;
    end
else
    for i = 1:lS
        dtmp = S(i).d;
        dtmp = circshift(dtmp,-(shift(i)-1));
        S(i).d = dtmp;
    end
end
