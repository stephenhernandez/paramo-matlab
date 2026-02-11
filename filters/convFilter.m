function df = convFilter(data,npts,zeroPhaseFlag)
%
% convFilter filter waveform data
%
% df = convFilter(data,npts,Fs,npoles,zeroPhaseFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Aug 29, 2023

%%
if nargin < 2
    npts = 3;
end

if nargin < 3
    zeroPhaseFlag = false;
end

%%
box = ones(npts,1)/npts;

% now filter
df = fftfilt(box,data);

% zero-phase filtering if requested

if zeroPhaseFlag
    df2 = flipud(data);
    df2 = fftfilt(box,df2);
    df2 = flipud(df2);
    df2 = cat(3,df,df2);
    df = mean(df2,3,"omitnan");
    df = squeeze(df);
end