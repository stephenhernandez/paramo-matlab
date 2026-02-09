function [nfft,d] = check_nfft(nfft,d)
if nargin < 2; d = []; end

%%
if isempty(d)
    %%
    if mod(nfft,2)
        %disp('removing a data point to make n even')
        nfft = nfft-1;
    end
    return;
end

%%
nfft = size(d,1);
if mod(nfft,2)
    %disp('removing a data point to make n even')
    d = d(1:nfft-1);
    nfft = nfft-1;
end