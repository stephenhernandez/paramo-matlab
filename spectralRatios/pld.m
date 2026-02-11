% projected landweber decon deconvolution test, sangay
clear; close all; clc;

newFs = 512; %seems like stuff goes weird above 400 sps, so that should be the max
secDur = 16;
maxIter = 100;

offsetDur1 = 0; %1.4;
offsetDur2 = offsetDur1 + 0.1;

cd ~/Desktop/;
try
    load('~/research/now/sangay/drumbeats/saga_drumbeat_data.mat');
catch
    try
        load('./igdata/saga_drumbeat_data.mat');
        catc
        fprintf(2,'couldnt load templates, abort');
        return;
    end
end

%
lfc = 0.25;
hfc = 4;
Sorig = loadWaveforms(datetime(2021,12,01),2,"SAGA",["HHZ"],"EC","",true,true); %,"~/research/now/sangay/drumbeats/");
S = intWaveforms(resampleWaveforms(filterWaveforms(detrendWaveforms(differentiateWaveforms(Sorig)),lfc,hfc),newFs));
Scut1 = cutWaveforms(S,t2-seconds(0),0,seconds(secDur));
amps2 = (pull(Scut1,'depmax') - pull(Scut1,'depmin'))*0.5;
[~,ampI] = sort(amps2,'descend');
t3 = t2(ampI);
amps3 = amps2(ampI);
Scut1 = Scut1(ampI);

S = intWaveforms(...
    resampleWaveforms(...
    filterWaveforms(...
    detrendWaveforms(...
    differentiateWaveforms(Sorig)),lfc,hfc),newFs));
Scut = cutWaveforms(S,t3-seconds(0),0,seconds(secDur));
Scut = detrendWaveforms(Scut);

d1 = pull(Scut1);
stack = pws(d1);
d0 = pull(Scut);
stack0 = pws(d0);
close all; figure(); plot(stack); zoom on;
figure(1); hold on; plot(stack0); zoom on;

[maxccp,delta] = doccFreqCircShift(taper(detrend(d0),0.5),true,[],[],newFs);

figure();
imagesc(squareform(maxccp)); zoom on; colorbar; caxis([-1 1]);

figure();
imagesc(squareform(delta)); zoom on; colorbar; caxis([-1000 1000]);

nG = size(d0,2);
Gshift = Gvdcc(nG);
W = (1/nG)*speye(nG);
raw_shifts = W*Gshift(1:end-1,:)'*delta;
figure(); plot(raw_shifts,'.'); zoom on;

raw_shifts = -round(raw_shifts);
shifted_data = apply_shifts(d0,raw_shifts);
stackShifted = pws(shifted_data);

figure(); plot(stackShifted); zoom on;

d2 = shifted_data(:,1:50);
plot_family(normalizeWaveforms(d2),(1:size(d2,2))',100,newFs); zoom on;

d2 = d0(:,1:50);
plot_family(normalizeWaveforms(d2),(1:size(d2,2))',100,newFs); zoom on;

Stest = Scut(1);
stackShifted = pws(detrend(shifted_data));
stackShifted = normalizeWaveforms(stackShifted);
Stest = dealHeader(Stest,stackShifted);
eGf = Stest; %cutWaveforms(S,t2(521)+seconds(offsetDur2),0,seconds(secDur));
egf = eGf(1);
egf = egf.d;
egf = detrend(egf);
Fs = 1./Stest(1).delta;
maxIter = 1000;
maxDur = 1/2;

we = abs(stackShifted);
we = we/sum(we);
amps4 = sum(we.*abs(shifted_data));

[prcntVector,peakCC,cc,eps,err_rms_prcnt] = cwiShifts(stackShifted,shifted_data,0,6,newFs,false,5,lfc,hfc);
figure(); semilogx(amps4(peakCC>=0.3),6*prcntVector(peakCC>=0.3)/100,'.'); zoom on; grid on;
figure(); imagesc(sign(shifted_data)'); colorbar; zoom on;

%%
loopN = size(shifted_data,2);
finalSTF = NaN(maxDur*newFs,loopN);
close all;
tic;
for i = 1:loopN%length(Stest)
    egf = eGf(1);
    egf = egf.d;
    egf = detrend((egf));
    %egf = circshift(egf,-round(0.25*newFs),1);
    egf = detrend(egf(round(0.25*newFs):end));
    legf = length(egf);
    tdum2 = (0:legf-1)';
    tdum = (0:legf-1)'./Fs;

    d = Scut(i).d;
    d = detrend((d(1:legf)));

    %     figure();
    %     ax(1) = subplot(311);
    %     plot(tdum,d,'LineWidth',2); grid on;
    %     ax(2) = subplot(312);
    %     plot(tdum,egf,'LineWidth',2); grid on;
    %     ax(3) = subplot(313);
    %     [cc,lags] = xcorr(d,egf);
    %     plot(lags/Fs,cc,'LineWidth',2); grid on;


    Gmax = 1./rssq(egf); %./max(abs(d));
    tau = 0.05/(Gmax); %rms(d)./rms(egf);
    stf = zeros(length(egf),maxIter);

    n = 1;
    ncols = 2;
    %fig = figure('units','normalized','position',[0 0 1 1]);
    %hold on;
    %ax1(1) = subplot(maxIter+1,ncols,n*ncols-1);
    %plot(stf(1:maxDur*newFs,1));
    %ax2(1) = subplot(maxIter+1,ncols,n*ncols);
    %synth = fftfilt((stf(1:maxDur*newFs,1)),egf);
    %plot(synth(1:length(d)));
    %hold on;
    %plot(d);

    mcc = NaN(maxIter,1);
    preverrormetric = 0;
    derror = 0;
    while n <= maxIter && derror <= 0
        n = n+1;
        thisStf = stf(:,n-1);
        synthOrig = filter((thisStf),1,egf);
        %synth = ifft(Synth,'symmetric');
        newStf = conv((tau*flipud(egf)),(d-synthOrig));
        newStf = fftshift(newStf);
        newStf = thisStf + newStf(1:length(thisStf));
        s_ = newStf(1:length(egf));
        %s_ = ifft(newStf,'symmetric');
        sI = tdum2 >= maxDur*newFs | s_ < 0;
        s_(sI) = 0;
        %s_ = s_/sum(s_);
        %newStf = s_;
        synth = filter(s_,1,egf);
        stf(:,n) = s_;

        %ax1(n) = subplot(maxIter+1,ncols,n*ncols-1);
        %plot(s_(1:maxDur*newFs,1));
        %ax2(n) = subplot(maxIter+1,ncols,n*ncols);
        %plot(synth(1:length(d)));
        %hold on;
        %plot(d);
        %mcc(n-1) = max(conv(d,flipud(synth)));
        %mcc(n-1) = rms(d - synth);
        thisMetric = rms(d - synth); %dot(d,synth)/norm(d)/norm(synth);
        derror = thisMetric - preverrormetric;
        if n == 2
            derror = 0;
        end
        mcc(n-1) = thisMetric;
        preverrormetric = thisMetric;
        disp(n);
    end
    %figure(1); hold on; plot(mcc,'.'); zoom on;
    %figure(2); hold on; plot(s_,'linewidth',2); zoom on;
    finalSTF(:,i) = s_(1:maxDur*newFs);
end
toc;

figure(); ax__(1) = subplot(211); plot(synth(1:length(d)),'linewidth',2);
hold on;
plot(d,'linewidth',2); zoom on; ax__(2) = subplot(212); plot(egf,'linewidth',2); zoom on; linkaxes(ax__,'x'); %plot(fftfilt(s_,egf),'o'); zoom on;

[prcntVector,peakCC,cc,eps,err_rms_prcnt] = cwiShifts(stackShifted,shifted_data,1,5,newFs,false,5,lfc,hfc);
we = abs(stackShifted);
we = we/sum(we);
amps4 = sum(we.*abs(shifted_data));
figure(); semilogx(amps4,prcntVector,'.'); zoom on; grid on;
figure(); semilogx(amps4(peakCC>=0.3),prcntVector(peakCC>=0.3),'.'); zoom on; grid on;

%%
figure(); plot(finalSTF); zoom on;
stf2 = finalSTF(:,2);
figure(); plot(stf2); zoom on;
for i = 1:size(finalSTF,2)
    stf_ = finalSTF(:,i);
    I = find(stf_==0,1);
    stf_ = stf_(I:end);
    finalSTF(1:length(stf_),i) = stf_;
end
stf2 = finalSTF(:,2);
figure(); plot(stf2); zoom on;
for i = 1:size(finalSTF,2)
    stf_ = finalSTF(:,i);
    I = find(stf_>0,1);
    if I > 1
        I = I -1;
    end
    stf_ = stf_(I:end);
    finalSTF(1:length(stf_),i) = stf_;
end
stf2 = finalSTF(:,2);
figure(); plot(stf2); zoom on;
for i = 1:size(finalSTF,2)
    stf_ = finalSTF(:,i);
    I = find(stf_==0,2);
    width(i) = I(end); maxSTF = max(stf_(1:width(i))); peak(i) = maxSTF;
end
figure(); plot(width,'.'); zoom on;
figure(); loglog(peak,width/newFs,'.'); zoom on;
