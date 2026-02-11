function [Hbu,Fbu] = freqOperator(npts,lfc,hfc,Fs,npoles)
if nargin < 5; npoles = 4; end
if nargin < 4;  Fs = 1; end

%%
Hd = zpkOperator(lfc,hfc,Fs,npoles);        %butterworth filter
[Hbu,Fbu] = freqz(Hd,npts,"whole",Fs);