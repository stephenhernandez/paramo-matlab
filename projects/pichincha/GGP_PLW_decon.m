cd ~/research/now/pichincha/;
clear; close all;
%cd ~/Desktop; 
load PichinchaRepeaters_4kEvents_14channels_2022_2023_displacement.mat;

%%
GfOrig = Gf;
diffFlag = true;
[Gf,maxAmp] = normalizeWaveforms(Gf,true,true,true); %amplitudes are modified!

nu = 1;
newFs = 256;
newLFC = 2;
newHFC = 16;

%
w = zeros(57358,1);
if diffFlag
    w = w(1:end-size(Gf,2),:);
end

linStack_ = w;
nTrace = size(Gf,1);
if diffFlag
    Gf = differentiateWaveforms(Gf);
end
Gf = filterWaveforms(taperWaveforms(detrendWaveforms(Gf),2*newFs),newLFC,newHFC); %amplitudes are modified!

for i = 1:nTrace
    Gf_ = Gf(i,:)';
    d_ = pull(Gf_);
    d_ = d_(:);
    dI = ~isfinite(d_);
    d_(dI) = 0;
    linStack_ = linStack_ + d_;
    w_ = exp(1j*(angle(hilbert(d_))));
    w = w + w_;
    disp(i);
end

w = abs(w/nTrace);
w = w.^nu;
linStack_ = linStack_/nTrace;

figure();
plot(w); zoom on;

figure();
plot(linStack_.*w); zoom on;

figure();
plot(linStack_); zoom on; hold on;
plot(linStack_.*w);

figure();
plot(w,'linewidth',1); grid on;
yyaxis right; plot(linStack_); zoom on;
hold on; plot(linStack_.*w,'k-','linewidth',2); grid on;

pwStack = linStack_.*w;
ampWeights = pwStack.^2;
ampWeights = ampWeights/sum(ampWeights);

allAmps = NaN(nTrace,1);
Gf = GfOrig;
if diffFlag
    Gf = differentiateWaveforms(Gf);
end
Gf = filterWaveforms(taperWaveforms(detrendWaveforms(Gf),2*newFs),newLFC,newHFC); %the amplitudes are fine!

for i = 1:nTrace
    Gf_ = Gf(i,:)';
    d_ = pull(Gf_);
    d_ = d_(:);
    dI = ~isfinite(d_);
    d_(dI) = 0;
    amp_ = sum(ampWeights.*(d_.^2));
    allAmps(i) = amp_;
    disp(i);
end

figure();
semilogy(pull(Gf(:,1),'ref'),allAmps,'.'); zoom on; grid on;

%%
close all;
nfft = 2^(nextpow2(Gf(1).npts+1));
Gf1 = Gf(:,1);

[sortAllAmps,sI] = sort(allAmps);
Gf1 = Gf1(sI);
d = pull(Gf1(1:3000)); %sortAllAmps<=2.5e5));
Gf1 = flipud(Gf1);

eGf = mean(normalizeWaveforms(detrend(d)),2,'omitnan');
egf_si = 2*newFs+150; egf_dur = 12*newFs; eGf2 = detrend(eGf(egf_si+20:egf_si+20+egf_dur-1));
E = fft(eGf2,nfft);
eGf2 = ifft(E,'symmetric');

figure();
plot(eGf2); zoom on;

tauInv = 1/rssq(conv(flipud(eGf2),eGf2,'full'));
Sall = zeros(nfft,nTrace);
tdum = (0:nfft-1)/newFs;
clear varRed;
varRed = NaN(nTrace,1);

si = zeros(nTrace,1);
mi = si;
ei = mi;
hd = ei;

wu = w(1:4097);
if diffFlag
    wu = wu(1:end-1);
end
wu = wu(2*newFs:end);

for i = 1:nTrace
    u = Gf1(i).d;
    u = u(egf_si:egf_si+egf_dur-1);
    %u = u.*wu; %experimental
    U = fft(u,nfft);
    u = ifft(U,'symmetric');


    S = fft(zeros(nfft,1),nfft);
    M = E.*S;
    R = U - M;

    S1 = ifft(S,'symmetric') + ifft(tauInv.*conj(E).*R,'symmetric');
    S1(S1<=0) = 0;
    S1(tdum >= 0.21) = 0;
    S2 = ifft(S,'symmetric')+S1;
    synth = conv(eGf2,S1,'full');
    synth = synth(1:nfft);

    % comment out
    % close all; figure(); plot(tdum,S1,'.-'); zoom on; grid on;
    % figure(); plot(tdum,synth); zoom on; grid on; hold on; plot(tdum,u); %plot(eGf2);
    % figure(); plot(tdum,[ifft(S,'symmetric') S1 S2],'.-'); zoom on;

    si_ = max([1 find(S2>0,1,'first')-1]);
    si(i,1) = si_;
    S2 = S2(si_:end);
    ei_ = 1+find(S2(2:end)==0,1,'first');
    ei(i,1) = ei_;
    S2 = S2(1:ei_);
    mi_ = find(S2==max(S2),1);
    mi(i,1) = mi_;

    S2(S2==0) = 1e-6;
    varRed(i) = sum(u.*synth)./(rssq(u).*rssq(synth));

    %disp(varRed(i));
    Sall(1:length(S2),i) = S2;
    hd_ = diff(interp1(cumsum(S2)/sum(S2),(0:length(S2)-1)',[0.25; 0.75]));
    hd(i,1) = hd_;

    % % comment out
    % S = fft(S1,nfft);
    % M = E.*S;
    % R = U - M;
    % S1 = ifft(S,'symmetric') + ifft(tauInv.*conj(E).*R,'symmetric');
    % S1(tdum <= si_/newFs) = 0;
    % S1(tdum >= 0.21) = 0;
    % S1(S1<0) = 0;
    % S2 = ifft(S,'symmetric')+S1;
    % synth = conv(eGf2,S2,'full');
    % synth = synth(1:nfft);
    % figure(); plot(tdum,[S1 S2],'.-'); zoom on; grid on;
    % figure(); plot(tdum,synth); zoom on; grid on; hold on; plot(tdum,u); %plot(eGf2);
    % figure(); plot(tdum,[ifft(S,'symmetric') S1 S2],'.-'); zoom on;
end

%
figure();
plot(varRed,'.'); zoom on; grid on;
ylabel('variance reduction');
ax = gca; ax.Color = [0.85 0.85 0.85];

[~,vrI] = sort(varRed,'descend');
figure();
plot((0:50-1)'/newFs,Sall(1:50,vrI(1:100)),'linewidth',1.5); zoom on; grid on;
ax = gca; ax.Color = [0.85 0.85 0.85];

tStart = pull(Gf1,'ref');
[tStart,tI] = sort(tStart);
timeToNext = seconds(diff(tStart));

figure();
semilogy(timeToNext,'.'); zoom on;

eiSort = ei(tI);
propMoment = sum(Sall)';
propMomentSort = propMoment(tI);

figure();
semilogy(tStart,propMomentSort,'.'); zoom on; grid on;
ax = gca; ax.Color = [0.85 0.85 0.85];

%
CCThresh = 0.4;
figure();
loglog(propMoment(varRed>=CCThresh),newFs./hd(varRed>=CCThresh),'.');
zoom on; title('moment-duration scaling is poor'); grid on;
ax = gca; ax.Color = [0.85 0.85 0.85];

tmpI = find(varRed>=CCThresh);
b_ = flipud(robustfit(log10(propMoment(tmpI)),log10(newFs./hd(tmpI))));
aSample = [min(log10(propMoment(tmpI))); ...
    log10(propMoment(tmpI(randsample(length(tmpI),100)))); ...
    max(log10(propMoment(tmpI)))];
yq = polyval(b_,aSample);
hold on; ll = loglog(10.^aSample,10.^yq,'.');
disp(b_);

%
figure();
subplot(121);
hdSeconds = hd/newFs;
loglog(sort(hdSeconds),1-((0:length(hdSeconds)-1)'/length(hdSeconds)),'.'); zoom on; grid on;
xlabel('half-duration'); title('power law');
ax = gca; ax.Color = [0.85 0.85 0.85];
subplot(122);
semilogy(sort(propMomentSort),1-((0:length(propMomentSort)-1)'/length(propMomentSort)),'.'); zoom on; grid on;
xlabel('moment ratio'); title('exponential');
ax = gca; ax.Color = [0.85 0.85 0.85];
