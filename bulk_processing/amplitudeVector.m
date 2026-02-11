function [ampVec,winlen] = amplitudeVector(d,Fs,dur,meanFlag,rmsFlag,exponent)
%
% amplitudeVector convert input data to rsam stream
%
% ampVec = amplitudeVector(d,Fs,dur,meanFlag,method,exponent)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
if nargin < 2
    Fs = 100;
end

if nargin < 3
    dur = 60;
end

if nargin < 4
    meanFlag = true;
end

if nargin < 5
    rmsFlag = false;
end

if nargin < 6
    exponent = 1;
end

if exponent ~= 1
    d = d.^exponent;
end

%%
if rmsFlag
    d = d.*conj(d);                 % amplitudes squared
else
    d = abs(d);
end

if Fs > 1
    Fs = round(Fs);
end

winlen = round(dur*Fs);
if meanFlag
    box = ones(winlen,1)/winlen;
    ampVec = fftfilt(box,d);        % sliding mean here (potentially slow)
else
    ampVec = medfiltSH(d,winlen);   % or med. filt here
end

if rmsFlag
    ampVec = sqrt(abs(ampVec));          % raiz cuadrada (square root of either mean or median square)
end
