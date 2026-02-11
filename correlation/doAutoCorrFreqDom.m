function [D,varargout] = doAutoCorrFreqDom(d,verboseFlag)
if nargin < 2
    verboseFlag = false;
end

%%
[winlen,n] = size(d);

%% pre-allocate data, get first ffts
nfft = nextpow2(winlen)+1;
nfft = 2^nfft;                  % nfft is always even
D = fft(d,nfft,1);              % convert to freq-dom

%%
if verboseFlag
    disp([winlen,n]);
end

%%
D = D.*conj(D);                 % perform auto-correlation
D = ifft(D,[],1,'symmetric');   % convert to time-domain
D = fftshift(D,1);
D = D(2:end,:);

%%
if nargout > 1
    m = nfft/2; %(size(D,1)+1)/2;
    lags = (-m+1:m-1)';
    varargout{1} = lags;
end