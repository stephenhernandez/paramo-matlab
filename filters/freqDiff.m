function dprime = freqDiff(d,Fs)
%
% freqDiff differentiate waveform in frequency domain
%
% dprime = freqDiff(d,Fs)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 1
    fprintf(2,'not enough input arguments, returning empty array\n');
    dprime = [];
    return;
end

if nargin < 2
    Fs = 1;
end

%%
rd = size(d,1);
nfft = 2^nextpow2(rd);

%%
D = fft(d,nfft);
iomega = diffOperator(nfft,Fs); % pad with zeros
D = D.*iomega;
dprime = ifft(D,'symmetric');   % make sure ifft is symmetric
dprime = dprime(1:rd,:);        % chop off previously padded zeros