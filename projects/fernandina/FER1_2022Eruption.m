clear; close all; clc;
warning off signal:findpeaks:largeMinPeakHeight

%
cd ~/Desktop/;
try
    load('~/research/now/fernandina/fer1_plantillas_hf.mat');
catch
    try
        load('./fer1_plantillas_hf.mat');
    catch
        fprintf(2,'couldnt load templates, abort');
        return;
    end
end

%%
maxN = 1e4;
tMaster = NaT(maxN,1);
pksMaster = NaN(maxN,1);
%threshLow = 0.12999;%2832157366
thresh = 0.2;
Fs = 100;
medthresh = 8;

%
dayStart = datetime(2020,01,12);
dayEnd = datetime(2020,01,12);

dayStart = datetime(2022,02,26);
dayEnd = datetime(2022,03,04);

dInc = 1;
days = (dayStart:dInc:dayEnd)';
ldays = length(days);

%
pwsFlag = 1;
if size(dz,2) > 1
    if pwsFlag
        dz = resample(pws(normalizeWaveforms(dz)),1,1);
        dn = resample(pws(normalizeWaveforms(dn)),1,1);
        de = resample(pws(normalizeWaveforms(de)),1,1);
    else
        dz = resample(mean(normalizeWaveforms(dz),2),1,1);
        dn = resample(mean(normalizeWaveforms(dn),2),1,1);
        de = resample(mean(normalizeWaveforms(de),2),1,1);
    end
end

kstnm = "FER1";
chanlist = ["BHZ";"BHN";"BHE"];

%%
templates = double(flipud(normalizeWaveforms([dz dn de])));
winlen = round(size(templates,1)/1)-1;

snippetDur = 10;
dcutmaxlength = snippetDur*Fs + 1;
daligned1 = zeros(dcutmaxlength,maxN);
daligned2 = zeros(dcutmaxlength,maxN);
daligned3 = zeros(dcutmaxlength,maxN);
templates = normalizeWaveforms(templates(end-winlen+1:end,:));

%%
lfc = 3;
hfc = 12;
npoles = 4;

%
si = 1;
nTot = 0;
for i = 1:ldays
    day_ = days(i);
    S = loadWaveforms(day_,dInc,kstnm,chanlist);
    if isnat(S(1).ref)
        fprintf(2,'no data for day: %s\n',datestr(day_));
        continue;
    end

    %
    Sf2 = intWaveforms(filterWaveforms(detrendWaveforms(differentiateWaveforms(S)),4,hfc,npoles));
    %Sf2 = transferWaveforms(S,3,12,4,Fs,'disp');
    Sf2 = resampleWaveforms(Sf2,Fs);
    Sf2 = syncWaveforms(Sf2);

    Sf = intWaveforms(filterWaveforms(detrendWaveforms(differentiateWaveforms(S)),lfc,hfc,npoles));
    Sf = resampleWaveforms(Sf,Fs);
    Sf = syncWaveforms(Sf);

    npts = pull(Sf,'npts');
    maxNpts = max(npts);


    %
    lS = length(Sf);
    ccnorm = zeros(maxNpts,lS);
    for k = 1:lS
        dlong = Sf(k).d;
        ld = npts(k);
        plantilla = templates(:,k);
        norms = sqrt(abs(fftfilt(ones(winlen,1),dlong.^2)));
        cc = fftfilt(plantilla,dlong);
        try
            %ccnorm(1:ld,k) = ccnorm(1:ld,k) + cc(1:ld)./norms(1:ld);
            ccnorm(1:ld,k) = cc(1:ld)./norms(1:ld);
        catch
            fprintf('something went wrong with length of files\n');
            continue;
        end
    end
    ccnorm = mean(ccnorm,2); %,'omitnan');

    %
    t = getTimeVec(Sf);

    %
    mad_ = mad(ccnorm,1);
    threshLow = medthresh*mad_;
    [pks,locs2] = findpeaks((ccnorm),'MINPEAKDISTANCE',winlen/4,'MINPEAKHEIGHT',threshLow);

    n_ = length(locs2);
    if ~n_
        fprintf(1,'no repeats found for day: %s\n',datestr(day_));
        continue;
    end

    %
    while locs2(1) < winlen
        if locs2(1) < winlen
            if n_ > 1
                locs2 = locs2(2:end);
                pks = pks(2:end);
                n_ = n_ - 1;
            else
                continue;
            end
        end
    end

    %
    t2 = t(locs2-winlen+1);
    tMaster(si:si+n_-1) = t2;
    pksMaster(si:si+n_-1) = pks;

    snippetShift = 1;
    for k = 1:lS
        Scut2 = cutWaveforms(Sf2(k),t2-seconds(snippetShift),0,seconds(snippetDur));
        daligned_ = double(pull(Scut2));

        %
        dlength = size(daligned_,1);
        if k == 1
            daligned1(1:dlength,si:si+n_-1) = daligned_;
        elseif k == 2
            daligned2(1:dlength,si:si+n_-1) = daligned_;
        else
            daligned3(1:dlength,si:si+n_-1) = daligned_;
        end
    end

    %
    fprintf(1,'found %d repeats for day: %s\n',n_,datestr(day_));
    nTot = nTot + n_;
    si = si + n_;
end

%
if ~nTot
    fprintf(2,'nothing processed\n');
    return;
end
    
%
tw = 40/winlen;
tMaster = tMaster(1:nTot);
pksMaster = pksMaster(1:nTot);
daligned1 = daligned1(:,1:nTot);
daligned2 = daligned2(:,1:nTot);
daligned3 = daligned3(:,1:nTot);

ampWindow = min([15*Fs dcutmaxlength]);
amps = 0.5*peak2peak(taper(detrend(daligned1(1:ampWindow,:)),tw))';
amps2 = 0.5*peak2peak(taper(detrend(daligned2(1:ampWindow,:)),tw))';
amps3 = 0.5*peak2peak(taper(detrend(daligned3(1:ampWindow,:)),tw))';
ampOrig = [amps amps2 amps3];

%%
tic;
close all;
dOrig = [daligned1; daligned2; daligned3];
ampRMS = 0.5*peak2peak(dOrig)'; %sqrt(2)*rms(abs(dOrig))'; %median([rms(daligned1(1:2000,:))' rms(daligned2(1:2000,:))' rms(daligned3(1:2000,:))'],2,'omitnan');
%goodI = (pksMaster >= 0.25 | ampRMS >= 8e2) & tMaster >= datetime(2021,12,01,12,00,00) & tMaster < datetime(2021,12,02,12,00,00);
goodI = pksMaster >= thresh; % & tMaster >= datetime(2021,12,01,12,00,00) & tMaster < datetime(2021,12,02,12,00,00);
%goodI = (tMaster >= datetime(2021,04,20) & tMaster < datetime(2021,05,01)) | (tMaster >= datetime(2021,12,01) & tMaster < datetime(2021,12,03));
%tmpStackAligned = plot_family(daligned1,1:size(daligned1,2),20,Fs);

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(tMaster(goodI),ampRMS(goodI),'o'); zoom on; grid on;
%figure(); imagesc((0:size(daligned1,1)-1)'/Fs,(1:size(daligned1,2))',sign(daligned1)'); zoom on; colorbar;

%d = normalizeWaveforms([normalizeWaveforms(daligned1); normalizeWaveforms(daligned2); normalizeWaveforms(daligned3)]);
dOrig = dOrig(:,goodI);
d = normalizeWaveforms(dOrig);
[pwsStack,~,lStack] = pws(d);
weights = abs(pwsStack)./sum(abs(pwsStack));
%weights = lStack/sum(abs(lStack));
amps2 = 0.8*2.2*sqrt(sum(weights.*(dOrig.^2))');
figure(1); hold on; semilogy(tMaster(goodI),amps2,'.'); zoom on; grid on;

figure('units','normalized','outerposition',[0 0 1 1]);
imagesc((0:size(d,1)-1)'/Fs,(1:size(d,2))',sign(d)'); zoom on; colorbar; colormap gray;

figure('units','normalized','outerposition',[0 0 1 1]);
ax1(1) = subplot(211);
tGood = tMaster(goodI);
semilogy(tGood(2:end),seconds(diff(tGood)),'.'); zoom on; grid on;
ax1(2) = subplot(212);
plot(tGood,(1:length(tGood))','.'); zoom on; grid on;
linkaxes(ax1,'x');

figure('units','normalized','outerposition',[0 0 1 1]);
plot(tGood,pksMaster(goodI),'.'); zoom on;

ratioWindow = 50;
if nTot < ratioWindow
    return;
end

%
[P,t_,~,stdT] = pratio(tGood,ratioWindow);
figure('units','normalized','outerposition',[0 0 1 1]);
ax2(1) = subplot(211);
semilogy(t_,P,'.'); zoom on;

ampGood = amps2; %medamps(goodI);
[PP,~,~,stdA] = pratio(amps2(2:end)./amps2(1:end-1),ratioWindow);
%PP = PP(2:end);
grid on;

ax2(2) = subplot(212);
semilogy(t_,PP,'.'); zoom on;
grid on;
linkaxes(ax2,'x');

figure('units','normalized','outerposition',[0 0 1 1]);
plot(t_,P./PP,'.'); zoom on;
toc;

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(tGood(2:end),ampGood(2:end)./ampGood(1:end-1),'.'); zoom on; grid on; title('variation in amplitude ratio from one event to the next');

figure('units','normalized','outerposition',[0 0 1 1]);
ll = scatter(seconds(diff(tGood)),ampGood(2:end)./ampGood(1:end-1),exp(log10(ampGood(2:end))),datenum(tGood(2:end)),'filled');
C = colorbar; zoom on; grid on; ax = gca; ax.XScale = 'log'; ax.YScale = 'log'; ll.MarkerFaceAlpha = 0.5; ll.MarkerEdgeColor = 'k'; ll.MarkerEdgeAlpha = 0.5;
C.TickLabels = datestr(C.Ticks);

tdiff1 = seconds(diff(tGood));
figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(tGood(3:end),(1*((tdiff1(2:end)./tdiff1(1:end-1)))),'.'); zoom on; grid on;


%%
%%
% queryDT = seconds(1);
% tquery1 = (dateshift(min(tGood),'end','minute'):queryDT:dateshift(max(tGood),'start','minute'))';
% smooth1 = interp1(datenum(tGood),(1:length(tGood))',datenum(tquery1));
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(tquery1(2:end),60*diff(smooth1),'.'); zoom on;
% figure(10); hold on; l2 = plot(tquery1(2:end),zpkFilter(60*diff(smooth1),-inf,1/600,1,1,1),'linewidth',2); l2.Color(4) = 0.75; zoom on; ax = gca; ax.YScale = 'log'; grid on;
% 
% tquery2 = (dateshift(min(tGood),'end','minute'):queryDT:dateshift(max(tGood),'start','minute'))';
% smooth2 = interp1(datenum(tGood),cumsum(ampGood),datenum(tquery2));
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(tquery2(2:end),0.1*diff(smooth2),'.'); zoom on;
% figure(11); hold on; l2 = plot(tquery2(2:end),zpkFilter(0.1*diff(smooth2),-inf,1/600,1,1,1),'linewidth',2); l2.Color(4) = 0.75; zoom on; ax = gca; ax.YScale = 'log'; grid on;
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(tquery2(2:end),zpkFilter(60*diff(smooth1),-inf,1/600,1,1,1).*zpkFilter(0.1*diff(smooth2),-inf,1/600,1,1,1),'linewidth',2); zoom on;
% ax = gca; ax.YScale = 'log'; grid on;
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(tquery2(2:end),cumsum(zpkFilter(60*diff(smooth1),-inf,1/600,1,1,1).*zpkFilter(0.1*diff(smooth2),-inf,1/600,1,1,1)),'linewidth',2); zoom on;
% ax = gca; ax.YScale = 'log'; grid on;
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(tGood,ampRMS(goodI)./amps2,'.'); zoom on;
% 
% %
% t = getTimeVec(Sf);
% figure('units','normalized','outerposition',[0 0 1 1]);
% wax(1) = subplot(411);
% lt = size(ccnorm(winlen:end),1);
% plot(t(1:lt),ccnorm(winlen:winlen+lt-1)); zoom on; hold on;
% plot(t(locs2-winlen+1),pks,'p');
% wax(2) = subplot(412);
% t = getTimeVec(Sf);
% plot(t,Sf2(1).d); zoom on; hold on; title('HHZ, 0.25 - 2 Hz.');
% wax(3) = subplot(413);
% t = getTimeVec(Sf(2));
% plot(t,Sf2(2).d); zoom on; hold on; title('HHN, 0.25 - 2 Hz.');
% wax(4) = subplot(414);
% t = getTimeVec(Sf(3));
% plot(t,Sf2(3).d); zoom on; hold on; title('HHE, 0.25 - 2 Hz.');
% linkaxes(wax,'x');
% 
% %
% tdum = seconds((0:winlen-1)/Fs);
% nsnippets = length(tMaster);
% snippetShiftN = snippetShift*Fs;
% Z = zeros(size(Sf(1).d));
% N = zeros(size(Sf(2).d));
% % dztemplate = mean(normalizeWaveforms(daligned1(1:winlen,:)),2); %/max(abs(daligned1(1:winlen)));
% % dntemplate = mean(normalizeWaveforms(daligned2(1:winlen,:)),2); %/max(abs(daligned1(1:winlen)));
% % detemplate = mean(normalizeWaveforms(daligned3(1:winlen,:)),2); %/max(abs(daligned1(1:winlen)));
% 
% dztemplate = mean(normalizeWaveforms(daligned1),2); %/max(abs(daligned1(1:winlen)));
% dntemplate = mean(normalizeWaveforms(daligned2),2); %/max(abs(daligned1(1:winlen)));
% detemplate = mean(normalizeWaveforms(daligned3),2); %/max(abs(daligned1(1:winlen)));
% 
% % dztemplate = pws(normalizeWaveforms(daligned1)); %/max(abs(daligned1(1:winlen)));
% % dntemplate = pws(normalizeWaveforms(daligned2)); %/max(abs(daligned1(1:winlen)));
% % detemplate = pws(normalizeWaveforms(daligned3)); %/max(abs(daligned1(1:winlen)));
% 
% % dztemplate = pws(normalizeWaveforms(daligned1(1:winlen,:))); %/max(abs(daligned1(1:winlen)));
% % dntemplate = pws(normalizeWaveforms(daligned2(1:winlen,:))); %/max(abs(daligned1(1:winlen)));
% % detemplate = pws(normalizeWaveforms(daligned3(1:winlen,:))); %/max(abs(daligned1(1:winlen)));
% 
% dztemplate = dztemplate/max(abs(dztemplate));
% dntemplate = dntemplate/max(abs(dntemplate));
% detemplate = detemplate/max(abs(detemplate));
% lZ = length(dztemplate);
% for i = 1:nsnippets
%     sign_ = 1;
%     if ccnorm(locs2(i)) < 0
%         sign_ = -1;
%     end
%     %ll = plot(wax(2),tMaster(i) + tdum,sign_*ampOrig(i,1)*dztemplate,'k','linewidth',2); ll.Color(4) = 0.75;
%     Z(locs2(i)-winlen+1-snippetShiftN:locs2(i)-winlen+1+lZ-1-snippetShiftN) = Z(locs2(i)-winlen+1-snippetShiftN:locs2(i)-winlen+1+lZ-1-snippetShiftN) + ampOrig(i,1)*dztemplate;
%     %ll = plot(wax(3),tMaster(i) + tdum,sign_*ampOrig(i,2)*dntemplate,'k','linewidth',2); ll.Color(4) = 0.75;
%     N(locs2(i)-winlen+1-snippetShiftN:locs2(i)-winlen+1+lZ-1-snippetShiftN) = N(locs2(i)-winlen+1-snippetShiftN:locs2(i)-winlen+1+lZ-1-snippetShiftN) + ampOrig(i,2)*dntemplate;
%     %ll = plot(wax(4),tMaster(i) + tdum,sign_*ampOrig(i,3)*detemplate,'k','linewidth',2); ll.Color(4) = 0.75;
% end
% %figure(15); xlim([datetime(2021,12,02,00,00,15) datetime(2021,12,02,00,00,45)]);
% 
% %
% figure('units','normalized','outerposition',[0 0 1 1]);
% t = getTimeVec(Sf2);
% axx(1) = subplot(211);
% plot(t,Sf2(1).d);
% grid on;
% hold on;
% plot(t,Z);
% axx(2) = subplot(212);
% plot(t,Sf2(1).d-Z);
% zoom on;
% grid on;
% linkaxes(axx,'xy');
% 
% %%
% %figure(); plot(t,ccnorm); zoom on;
% %clear S* cc* norms t dlong daligned_
% %save('drumbeat_search_results_thresh20_hf');