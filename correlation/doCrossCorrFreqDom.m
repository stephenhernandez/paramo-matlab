function [data1,varargout] = doCrossCorrFreqDom(data1,data2,verboseFlag)
%
% CC = doCrossCorrFreqDom(data1,data2,verboseFlag)
% [CC,lags] = doCrossCorrFreqDom(data1,data2,verboseFlag)
%

if nargin < 3
    verboseFlag = false;
end

%%
[winlen1,n1] = size(data1);
[winlen2,n2] = size(data2);

if winlen1 ~= winlen2
    if verboseFlag
        fprintf("durations between two input matrices are not the same.\n");
        fprintf("using the larger of the two.\n");
        fprintf("\n");
    end
end

n = n1;
if n1 ~= n2
    if verboseFlag
        fprintf("number of traces between two input matrices are not the same.\n");
        fprintf("using the smaller of the two.\n");
        fprintf("\n");
    end

    if n1 > n2
        % data1 too long, truncate to save of data2
        data1 = data1(:,1:n2);
        n = n2;
    elseif n2 > n1
        % data2 too long, truncate to save of data1
        data2 = data2(:,1:n1);
    end
end

%% pre-allocate data, get first ffts
winlen = max([winlen1 winlen2]); %size(data1);
nfft = nextpow2(winlen)+1;
nfft = 2^nfft;                                      %nfft is always even
data1 = fft(data1,nfft,1);
data2 = conj(fft(data2,nfft,1));

%%
if verboseFlag
    fprintf("%d %d\n",winlen,n);
end

%%
data1 = data1.*data2;
data1 = ifft(data1,[],1,"symmetric");
data1 = fftshift(data1,1);
data1 = data1(2:end,:);

%%
if nargout > 1
    m = (size(data1,1)+1)/2;
    lags = (-m+1:m-1)';
    varargout{1} = lags;
end