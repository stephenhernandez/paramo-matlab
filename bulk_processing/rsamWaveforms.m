function [S,ampVecOrig,dOrig,winlen,iStart] = rsamWaveforms(S,varargin)
%
% rsamWaveforms return in structure rsam of original waveforms
%
% S = rsamWaveforms(S,dur,lfc,hfc,npoles,meanFlag,method,diffFlag,zeroPhaseFlag,taperWidth)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
nVarargin = length(varargin);
functionDefaults = {...
    60,...      % dur
    -inf,...    % lfc
    -inf,...    % hfc
    4,...       % npoles
    true,...    % meanFlag
    true,...    % rmsFlag
    false,...   % diffFlag
    false,...   % zeroPhaseFlag
    false,...   % tw
    false,...   % medfiltFlag
    1};         % exponent (default 1)

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[dur,lfc,hfc,npoles,meanFlag,rmsFlag,diffFlag,zeroPhaseFlag,tw,medfiltFlag,exponent] = deal(optsToUse{:});

%%
if isempty(dur)
    dur = 60;
end

if isempty(lfc)
    lfc = -inf;
end

if isempty(hfc)
    hfc = -inf;
end

if isempty(npoles)
    npoles = 4;
end

if isempty(meanFlag)
    meanFlag = true;        % sliding mean (default) or sliding median?
end

if isempty(rmsFlag)
    rmsFlag = false;
end

if isempty(diffFlag)
    diffFlag = false;
end

if isempty(zeroPhaseFlag)
    zeroPhaseFlag = false;
end

if isempty(tw)
    tw = false;
end

if isempty(medfiltFlag)
    medfiltFlag = false;
end

if isempty(exponent)
    exponent = 1;
end

%%
if medfiltFlag
    S = medfiltWaveforms(S,medfiltFlag,true);
    S = demeanWaveforms(S);
end

if diffFlag
    S = differentiateWaveforms(S);
end

if isfinite(lfc)
    S = detrendWaveforms(S);
end

if tw
    S = taperWaveforms(S,tw);
end

cornersfin = isfinite([lfc hfc]);
if any(cornersfin)
    S = filterWaveforms(S,lfc,hfc,npoles,tw,zeroPhaseFlag);
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
if lS > 1
    S = syncWaveforms(S);
end

for i = 1:lS
    S_ = S(i);
    dOrig = S_.d;

    ref = S_.ref;
    delta = S_.delta;
    Fs = 1/delta;
    if Fs >= 1
        Fs = round(Fs);
    end

    [ampVecOrig,winlen] = amplitudeVector(dOrig,Fs,dur,meanFlag,rmsFlag,exponent);
    tref = dateshift(ref,"end","minute");
    iStart = t2i(tref,ref,delta);
    %ampVec = downsample(ampVecOrig,winlen,iStart-1);
    ampVec = ampVecOrig(iStart:winlen:end);
    S_ = dealHeader(S_,ampVec,1/dur,tref);
    S(i) = S_;
end
S = reshape(S,sizeS);
