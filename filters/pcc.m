%function pcc
clear; close all; clc;
S = loadWaveforms(datetime(2014,05,03),2,"PINO","SHZ","EC","",0,1);
Sf = filterWaveforms(differentiateWaveforms(S),0.75,12);
Scut = cutWaveforms(Sf,datetime(2014,05,03,20,33,06),0,10);
Scut = detrendWaveforms(Scut);

%%
master = Sf.d;
analytic1 = hilbert(master);
phi = angle(analytic1);
masterSin = abs(sin(phi));
masterCos = abs(cos(phi));

%%
f = Scut.d;
analytic = hilbert(Scut.d);
h = imag(analytic);

%%
phase = angle(analytic);
%phi = phase;
psi = phase;
templateSin = abs(sin(psi));
templateCos = abs(cos(psi));

%%
tic;
winlen = length(psi);
cossum = conv(masterCos,templateCos,'valid')/winlen; toc;
sinsum = conv(masterSin,templateSin,'valid')/winlen; toc;
lhs = cossum + sinsum; toc;


%%
% tic;
% phisum = (phisum - sumpsi)/2;
% toc;
% 
% %%
% tic;
% cpcc = (abs(cos(phisum)) - abs(sin(phisum)));
% toc;

%%
% tic;
% N = length(analytic);
% sum(abs(cos((phi - psi)./2)) - abs(sin((phi - psi)./2)))/N;
% toc;

%%


%%
% %w = pcc(u,v)
% m = size(u,1);
% n = size(v,1);
% 
% wn = m + n -1;
% w = NaN(wn,1);
% for k = 1:wn
%     %disp(k)
%     %w_ = 0;
%     %for j = max(1,k+1-n):min(k,m)
%     j = max(1,k+1-n):min(k,m);
%     %w_ = sum(u(j).*v(k-j+1));
%     %    gam = u(j)-v(k-j+1);
%     %    w_ = w_ + (abs(cos(gam/2)) - abs(sin(gam/2)));
%     %end
%     w(k) = sum(u(j).*v(k-j+1)); %w_;
% end