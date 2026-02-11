function [S,p] = dcDetrendWaveforms(S)
%
% discontinuous detrendWaveforms returns structure with detrended waveforms
%
% [S,p] = detrendWaveforms(S)
%
% S: structure with detrended waveforms
% p: polynomial coefficients used to model linear trend
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%% NOTE!!!
% this code hasnt been changed yet, its it need to be modified
% do _not_ use
%%

%%
S = S(:);
lS = size(S,1);
for i = 1:lS
    if ~isnat(S(i).ref)
        if S(i).gapFlag
            %%
            npts = S(i).npts;
            d = S(i).d;
            gapInfo = S(i).gapInfo;
            
            %%
            nGaps = size(gapInfo,1) + 1;
            goodStarts = [1; sum(gapInfo,2)];
            goodEnds = [gapInfo(:,1)-1; npts];
            for kk = 1:nGaps
                goodStart = goodStarts(kk);
                goodEnd = goodEnds(kk);
                d_2 = d(goodStart:goodEnd);
                d(goodStart:goodEnd) = detrend(demean(d_2));
            end
            
            %%
            [minVals,maxVals,meanVals] = minmaxmean(d);
            S(i).depmin = minVals;
            S(i).depmax = maxVals;
            S(i).depmen = meanVals;
        else
            %%
            d = S(i).d;
            npts = S(i).npts;
            d = demean(d);
            t = (0:npts-1)';
            p = polyfit(t,d,1);
            dSynth = p(1)*t + p(2);
            d = d - dSynth;
            d = demean(d);
            S(i).d = d;
            
            %%
            [minVals,maxVals,meanVals] = minmaxmean(d);
            S(i).depmin = minVals;
            S(i).depmax = maxVals;
            S(i).depmen = meanVals;
        end
    end
end