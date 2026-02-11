function S = clipWaveforms(S,clipLevel,rmsFlag)
%
% clipWaveforms return structure with clipped waveforms
%
% S = clipWaveforms(S,clipLevel,rmsFlag)
% clipLevel: values above and below clipLevel are set to clipLevel
% rmsFlag: if rmsFlag true, clipLevel is set to clipLevel*rmsValue
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    clipLevel = 1e3;
end

if nargin < 3
    rmsFlag = false;
end

%%
lS = size(S,1);
if rmsFlag
    % if rmsFlag true, cliplevel is number of rms multiples to clip at
    for i = 1:lS
        if ~isnat(S(i).ref)
            d = S(i).d;
            gI = isfinite(d);
            d_ = d(gI);
            rmsValue = rms(d_);
            sgn = sign(d_);
            dI = abs(d_) >= clipLevel*rmsValue;
            d_(dI) = rmsValue*clipLevel.*sgn(dI);
            d(gI) = d_;
            S(i).d = d;
            
            %%
            [minVals,maxVals,meanVals] = minmaxmean(d);
            S(i).depmin = minVals;
            S(i).depmax = maxVals;
            S(i).depmen = meanVals;
        end
    end
else
    for i = 1:lS
        if ~isnat(S(i).ref)
            d = S(i).d;
            sgn = sign(d);
            dI = abs(d) >= clipLevel;
            d(dI) = clipLevel.*sgn(dI);
            S(i).d = d;
            
            %%
            [minVals,maxVals,meanVals] = minmaxmean(d);
            S(i).depmin = minVals;
            S(i).depmax = maxVals;
            S(i).depmen = meanVals;
        end
    end
end
