function [H,f] = cmplxResp(nfft,zeros,poles,constant,Fs)
% function [H,f] = cmplxResp(x,zero,pole,constant,Fs)
%
% transfer function frequency-domain response, given poles, zeros, and gain
% nfft must be even, and Fs is your sampling interval
%
% want the inverse transform for freq.-dom. deconvolution and such?
% try: [Hinv,f] = cmplxResp(z,pole,zero,1/constant,Fs)
%
% this code initially started as a transliterated version of getran.c in
% sac101.6a, where getran is called by polezero, which is called by dseis,
% which is called by transfer.c, which is parsed by xtransfer.
% basically, lots of complex algebra.
%
% but then matlab (via `freqs') taught me a faster way to accomplish the
% task using polyval and poly (see code).
%
% a much (much!) slower option (not utlized here) would be to loop:
% H(s) = H(1j*2*pi*f) = prod(s-zero)/prod(s-pole)
% with `zero' and `pole' as above
% perhaps seeing it written in this form helps to clarify the algorithm

%nfft = check_nfft(nfft);
H = complex(ones(nfft,1));

%% # positive frequencies
oddFlag = mod(nfft,2);
if oddFlag
    nfreq = (nfft+1)/2;
else
    nfreq = 1 + (nfft/2);
end

%%
f = (0:nfft-1)'*Fs/nfft;    % df = Fs/nftt
s = 1j*2*pi*f(1:nfreq);     % Laplace complex variable

% first half (Horner's Method) what happens if poly(zero) or poly(pole)
% returns complex data?
H(1:nfreq) = polyval(poly(zeros),s)./polyval(poly(poles),s);

% exploit symmetry
if oddFlag
    H(nfreq+1:end) = flipud(conj(H(2:nfreq)));
else
    H(nfreq+1:end) = flipud(conj(H(2:nfreq-1)));
end

% the DC should be 0 (no response at zero frequency)
H(1) = 0;

if ~oddFlag
    % the Nyquist should be real
    H(nfreq) = abs(H(nfreq));
end

% dunzo
H = constant*H;


% partial implementation (for pedagogic purposes) of old code that i pulled
% from SAC's getran.c
% for i = 2:nfreq
%     trn = 1.0;
%     tin = 0.0;
%     for j = 1:length(zero)
%         tr = -real(zero(j));
%         ti = 2*pi*f(i) - imag(zero(j));
%         tr0 = trn*tr - tin*ti;
%         ti0 = trn*ti + tin*tr;
%         trn = tr0;
%         tin = ti0;
%     end
%
%     trd = 1.0;
%     tid = 0.0;
%     for j = 1:length(pole)
%         tr = -real(pole(j));
%         ti = 2*pi*f(i) - imag(pole(j));
%         tr0 = trd*tr - tid*ti;
%         ti0 = trd*ti + tid*tr;
%         trd = tr0;
%         tid = ti0;
%     end
%
%     fac = constant / (trd^2 + tid^2);
%     H(i) = fac*complex(trn*trd + tin*tid,trd*tin - trn*tid);
% end