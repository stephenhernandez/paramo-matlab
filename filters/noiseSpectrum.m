clearvars -except zz;
close all; 
%["BRUN";"BULB";"BPAT";"BMAS";"BBIL"];
%S = loadWaveforms(datetime(2016,11,01),2,"BBIL",["BHZ";"BHN";"BHE"]);
S = loadWaveforms(datetime(2010,04,10),1,"BRUN",["BHZ";"BHN";"BHE"]);
S = loadWaveforms(datetime(2022,10,01),1,"BREF",["BHZ";"BHN";"BHE"]);
%S = loadWaveforms(datetime(2015,04,10),2,"PINO",["SHZ";"SHN";"SHE"]);

%%
Fs = 50;
secDur = 120;
lowestFreq = 0.01;
highestFreq = 25;
nfreqs = 1+1024*4;
freqVec = logspace(log10(lowestFreq),log10(highestFreq),nfreqs)';

winlen = secDur.*Fs;
nfft = secDur*Fs;
interpMethod = "pchip";

nOverlap = 0.5;
nfft3 = winlen;
nfft4 = (nfft3/2)+1;
tbw = 3;
[e,v] = dpss(nfft,tbw,10);
fxx = (0:Fs/nfft3:Fs/2)';

%%
S = differentiateWaveforms(S);
S = resampleWaveforms(S,Fs);
S = interpolateWaveforms(S);
S = syncWaveforms(S);
%S = scaleWaveforms(S,1/4);

minB = min(pull(S,'ref'));
maxE = max(pull(S,'ref')+pull(S,'e'));

tStart_ = dateshift(minB,'end','second');
tEnd_ = dateshift(maxE,'start','second');

dur_ = seconds(tEnd_-tStart_);
S = cutWaveforms(S,tStart_,0,dur_);

refs = pull(S,'ref');
badI = isnat(refs);
S(badI) = [];
lS = length(S);

%% experimental time-domain normalization
% d = double(pull(S));
% normers1 = max(abs(d),[],2,"omitnan");
% d = d./normers1;
% S(1).d = (d(:,1));
% S(2).d = (d(:,2));
% S(3).d = (d(:,3));

%%
dcut = [];
for j = 1:lS
    d_ = double(pull(S(j)));
    dcut_ = cutWindows(d_,winlen,nOverlap,true);
    dcut = cat(3,dcut,dcut_);
end

%%
nWindows = size(dcut,2);
nzz = NaN(nfft3,nWindows);
pzz = nzz;
pnz = nzz;
pze = nzz;
pne = nzz;

w = kaiser(winlen,10);

%%
for i = 1:nWindows
    disp(i);
    dcut_ = squeeze(dcut(:,i,:));
    %dcut_ = dcut_./median(rssq(dcut_));
    
    normers1 = max(abs(dcut_),[],2);
    dW = 0;
    timeDomain = false;
    verboseFlag = false;

    dcut_ = dcut_.*w;
    
    %dcut_ = dcut_./normers1;
    
    %dcut_ = fft(dcut_,winlen);
%     freqNormers = max(dcut_,[],2,"omitnan");
%     dcut_ = ifft(dcut_,'symmetric');
%     freqSmooth = zpkFilter(abs(freqNormers),-inf,1/15,1,1,1);

    dcut_ = [dcut_ dcut_(:,2)];
    dcut_2 = [repmat(dcut_(:,1),1,2) repmat(dcut_(:,3),1,2)];
    dcut_(:,3) = dcut_(:,1);

    D = fft(dcut_,winlen);
    D2 = fft(dcut_2,winlen);
    
%     freqNormers = max(D,[],2,"omitnan");
%     freqSmooth = zpkFilter(abs(freqNormers),-inf,1/15,1,1,1);
    %D(:,2:end) = D(:,2:end)./freqSmooth;
    %D2(:,2:end) = D2(:,2:end)./freqSmooth;

%     D_ = D(:,2:end);
%     D_ = ifft(D_,'symmetric');
%     D_ = D_./normers1;
%     D_ = fft(D_,winlen);
%     %D_ = fdWhiten(D_,-inf,-inf,dW,Fs,timeDomain,verboseFlag);
%     D(:,2:end) = D_;
% 
%     D_ = D2(:,2:end);
%     D_ = ifft(D_,'symmetric');
%     D_ = D_./normers1;
%     D_ = fft(D_,winlen);
%     %D_ = fdWhiten(D_,-inf,-inf,dW,Fs,timeDomain,verboseFlag);
%     D2(:,2:end) = D_;

%     D = fft(dcut_,winlen);
%     freqNormers = max(D,[],2,"omitnan");
%     freqSmooth = zpkFilter(abs(freqNormers),-inf,1/15,1,1,1);
%     D = D./freqSmooth;

    D = D.*conj(D2);

    %D = fftshift(ifft(D.*conj(D2)));
    %[P,f] = pmtm(dcut_,e,v,winlen,Fs);
    %[P,f] = pwelch(dcut_(:,1),w,0,winlen,Fs);
    %PZN = sqrt(P(:,1).*P(:,2)); %cpsd(dcut_(:,1),dcut_(:,2),w,0,winlen,Fs);
    %PZE = sqrt(P(:,1).*P(:,3)); %cpsd(dcut_(:,1),dcut_(:,3),w,0,winlen,Fs);
    %sPNE = sqrt(P(:,2).*P(:,3)); %cpsd(dcut_(:,2),dcut_(:,3),w,0,winlen,Fs);

    %[P,f] = pwelch(dcut_,w,nfft/4,winlen,Fs);

    %Pzz = P(:,1);
    %Pzz = cpsd(dcut_(:,1),dcut_(:,1),w,nfft/4,winlen,Fs);
    %ZN = cpsd(dcut_(:,1),dcut_(:,2),w,nfft/4,winlen,Fs);
    %ZE = cpsd(dcut_(:,1),dcut_(:,3),w,nfft/4,winlen,Fs);
    %NE = cpsd(dcut_(:,2),dcut_(:,3),w,nfft/4,winlen,Fs);

    %Pzz = P(:,1);
    %Pzn = sqrt(P(:,1).*P(:,2));
    %Pze = sqrt(P(:,1).*P(:,3));
    %Pne = sqrt(P(:,2).*P(:,3));
    %nzz_ = Pzz - (abs(ZN).*abs(ZE)./abs(NE));

    D = D./Fs./nfft3;
    D = abs(D);
    D = zpkFilter(D,-inf,1/40,1,1,1);

    %hold on; semilogy(abs(D.*exp(1j*phase)),'color',[0.25 0.25 0.25]); zoom on; grid on;
    %D = D(1:nfft4,:);
    %D(2:end-1,:) = 2*D(2:end-1,:);

    %nzz_ = (D(:,1)) - (D(:,2).*(D(:,3)./D(:,4)));
    %nzz_ = P(:,1) - abs(P(:,2)).*abs(P(:,3))./abs(P(:,4));

    %%
    %    nzz(:,i) = nzz_;
    %     pzz(:,i) = P(:,1);
    %     pzn(:,i) = abs(PZN);
    %     pze(:,i) = abs(PZE);
    %     pne(:,i) = abs(PNE);

    pzz(:,i) = D(:,1);
    pnz(:,i) = D(:,2);
    pze(:,i) = D(:,3);
    pne(:,i) = D(:,4);
end

%%
pureSignalHaHa = pnz.*pze./pne;

pzz2 = pzz(1:nfft4,:);
pnz2 = pnz(1:nfft4,:);
pze2 = pze(1:nfft4,:);
pne2 = pne(1:nfft4,:);
pureSignalHaHa = pureSignalHaHa(1:nfft4,:);

pzz2(2:end-1,:) = 2*pzz2(2:end-1,:); %D = zpkFilter(D,-inf,1/10,1,4,1);
pnz2(2:end-1,:) = 2*pnz2(2:end-1,:);
pze2(2:end-1,:) = 2*pze2(2:end-1,:);
pne2(2:end-1,:) = 2*pne2(2:end-1,:);
pureSignalHaHa(2:end-1,:) = 2*pureSignalHaHa(2:end-1,:);

nzz = pzz2 - pureSignalHaHa;
badI = abs(pzz2) < abs(pureSignalHaHa);

% nzz = sqrt(Fs.*nfft3.*nzz);
% pzz2 = sqrt(Fs.*nfft3*pzz2);

% pnz2 = sqrt(Fs.*nfft3*pnz2);
% pze2 = sqrt(Fs.*nfft3*pze2);
% pne2 = sqrt(Fs.*nfft3*pne2);
% pureSignalHaHa = sqrt(Fs.*nfft3.*pureSignalHaHa);

nzz = sqrt(Fs.*nfft3.*abs(nzz)); %nzz = zpkFilter(nzz,-inf,1/10,1,4,1);
pzz2 = sqrt(Fs.*nfft3*abs(pzz2)); %pzz2 = zpkFilter(pzz2,-inf,1/10,1,4,1);
pnz2 = abs(pnz2); %pnz2 = zpkFilter(pnz2,-inf,1/10,1,4,1);
pze2 = abs(pze2); %pze2 = zpkFilter(pze2,-inf,1/10,1,4,1);
pne2 = abs(pne2); %pne2 = zpkFilter(pne2,-inf,1/10,1,4,1);
pureSignalHaHa = sqrt(Fs.*nfft3.*abs(pureSignalHaHa)); %pureSignalHaHa = zpkFilter(pureSignalHaHa,-inf,1/10,1,4,1);

% pnz2 = sqrt(Fs.*nfft3*abs(pnz2));
% pze2 = sqrt(Fs.*nfft3*abs(pze2));
% pne2 = sqrt(Fs.*nfft3*abs(pne2));
% nzz = abs(nzz);
% pzz2 = abs(pzz2);
% pureSignalHaHa = abs(pureSignalHaHa);

nzz(badI) = NaN;
pureSignalHaHa(badI) = NaN;

%%
close all; 
figure(); 
loglog(fxx,median(pzz2,2,"omitnan"),'linewidth',4); hold on; grid on; zoom on; 
loglog(fxx,median(pnz2,2,"omitnan"),'linewidth',4); 
loglog(fxx,median(pze2,2,"omitnan"),'linewidth',4); 
loglog(fxx,median(pne2,2,"omitnan"),'linewidth',4);
loglog(fxx,(median(pnz2,2,"omitnan").*median(pze2,2,"omitnan")./median(pne2,2,"omitnan")),'k','linewidth',4); 
loglog(fxx,median(pureSignalHaHa,2,"omitnan"),'color',[0.5 0.5 0.5],'linewidth',4); 

legend(["Pzz";"Pnz";"Pze";"Pne";"Pure Signal 1";"Pure Signal Orig"]);

figure(); 
loglog(fxx,median(pnz2,2,"omitnan").*median(pze2,2,"omitnan"),'linewidth',4); hold on; zoom on; grid on; 
loglog(fxx,median(pne2,2,"omitnan"),'linewidth',4); 
loglog(fxx,median(pnz2,2,"omitnan").*median(pze2,2,"omitnan")./median(pne2,2,"omitnan"),'k','linewidth',4);

n2 = median(pzz2,2,"omitnan") - (median(pnz2,2,"omitnan").*median(pze2,2,"omitnan")./median(pne2,2,"omitnan"));

figure(); 
loglog(fxx,n2,'linewidth',4); zoom on; grid on; hold on; 
loglog(fxx,median(pzz2,2,"omitnan"),'linewidth',4);
loglog(fxx,median(nzz,2,"omitnan"),'linewidth',4);
loglog(fxx,(median(pnz2,2,"omitnan").*median(pze2,2,"omitnan")./median(pne2,2,"omitnan")),'k','linewidth',4); 
loglog(fxx,median(pureSignalHaHa,2,"omitnan"),'color',[0.5 0.5 0.5],'linewidth',4); 
%figure(); loglog(fxx,-1 + (median(pzz2,2,"omitnan")./n2),'linewidth',4); zoom on; grid on; %hold on; loglog(fxx,median(pzz2,2,"omitnan"),'linewidth',4);
