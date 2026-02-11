function df = zpkFilter(data,lfc,hfc,Fs,npoles,zeroPhaseFlag)
%
% zpkFilter filter waveform data
%
% df = zpkFilter(data,lfc,hfc,Fs,npoles,zeroPhaseFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

% modified 12 MAY 2023, different algorithm to implement zero phase

%%
df = data;
if nargin < 2
    fprintf('need at least two inputs, doing nothing\n');
    return;
end

if nargin < 3
    hfc = -inf;
end

if nargin < 4
    Fs = 100;
end

if nargin < 5
    npoles = 4;
end

if nargin < 6
    zeroPhaseFlag = false;
end

%%
cornersfin = isfinite([lfc hfc]);
if ~any(cornersfin)
    fprintf(2,'no valid corners input, doing nothing\n')
    return;
end

%if high pass filtering, remove trend...
if cornersfin(1)
    df = detrend(df);
end

if zeroPhaseFlag
    warning off signal:filtfilt:ParseSOS
    [~,sos,k] = zpkOperator(lfc,hfc,Fs,npoles);
    df = filtfilt(sos,k,df);
else
    Hd = zpkOperator(lfc,hfc,Fs,npoles);
    df = filter(Hd,df);
end
