function [Hd,sos,k] = zpkOperator(lfc,hfc,Fs,npoles)
%
% zpkOperator generate butterworth filter object
%
% Hd = zpkOperator(lfc,hfc,Fs,npoles)
%
% lfc: low-frequency corner
% hfc: high-frequency corner
% Fs: sampling rate
% npoles: number of poles
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    hfc = -inf;
end

if nargin < 3
    Fs = 1;
end

if nargin < 4
    npoles = 4;
end

%% get poles and zeros and convert to second-order section filter object
Fny = 0.5*Fs;       % Nyquist frequency
lfcn = lfc/Fny;
hfcn = hfc/Fny;

%% potentially dangerous, may give wrong impression about filter applied
if hfcn > 1
    fprintf('hfcn is greater than 1, ignoring value, assuming high-pass until nyquist...\n');
    hfcn = -inf;
    hfc = hfcn;
end

%%
cornersfin = isfinite([lfc hfc]);
if ~any(cornersfin)
    error('error, no valid corners input\n');
end

%%
if lfc < 0 % lfc is negative, assumes that hfc is positive (but butter will return error if something goes wrong)
    % low-pass: f < hfc
    [z,p,k]	= butter(npoles,hfcn);
    [sos,k] = zp2sos(z,p,k);
    Hd = dfilt.df2sos(sos,k);
    return;
end

if hfc > lfc % since lfc is positive => hfc is positive
    % band-pass, lfc<f<hfc
    [z,p,k]	= butter(npoles,[lfcn hfcn]);   
    [sos,k] = zp2sos(z,p,k);
    Hd = dfilt.df2sos(sos,k);
    return;
end

if hfc < 0 % hfc is negative, but lfc is positive
    % high-pass: f > lfc
    [z,p,k]	= butter(npoles,lfcn,'high');
    [sos,k] = zp2sos(z,p,k);
    Hd = dfilt.df2sos(sos,k);
    return;
end

if lfc > hfc % neither lfc nor hfc are negative, but hfc is less than lfc
    % band-stop hfc < f < lfc
    [z,p,k]	= butter(npoles,[hfcn lfcn],'stop');
    [sos,k] = zp2sos(z,p,k);
    Hd = dfilt.df2sos(sos,k);
    return;
end


