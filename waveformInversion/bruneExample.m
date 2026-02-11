clear; close all; clc;
pmtmFlag = true;
diffFlag = true;
boatwright = true;
curveFitFlag = false;
lfc = 2.;
hfc = 8;
kstnm = "CHL2";

S = loadWaveforms(datetime(2014,10,20),1,kstnm,"HHZ");
if diffFlag
    S = differentiateWaveforms(S);
end
S = filterWaveforms(S,lfc,hfc);
d = S.d;
t = getTimeVec(S);
tw = 0.1;

%tI = t >= datetime([2014 10 28 01 39 12.25]) & t < datetime([2014 10 28 01 39 14.75]);
%tI = t >= datetime([2014 10 19 12 27 47]) & t < datetime([2014 10 19 12 27 51]);
%tI = t >= datetime([2014 10 20 08 41 50]) & t < datetime([2014 10 20 08 41 55]);
%tI = t >= datetime([2014 10 20 08 41 36]) & t < datetime([2014 10 20 08 41 42]);
tI = t >= datetime([2014 10 20 08 42 56]) & t < datetime([2014 10 20 08 43 00]);

winlen = sum(tI);
eGf = double(d(tI));

%%
tic;
disp('getting data');
S = loadWaveforms(datetime(2014,10,10),14,kstnm,"HHZ");
if diffFlag
    S = differentiateWaveforms(S);
end
dorig = detrend(double(S.d));
S = filterWaveforms(S,lfc,hfc);
d = double(S.d);
toc;

%%
Fs = 1/S(1).delta;
t = getTimeVec(S);
clear S

%%
%perform normalized cross correlation
disp('perform normalized cross correlation')
tic;
box = d.^2;
box = conv(ones(winlen,1),box);
box = sqrt(box);
badI = box < 1 | ~isfinite(box);
box(badI) = 1;
toc;

tic;
disp('getting cross correlation');
cc = conv(flipud(eGf),d);
toc;

normcc = cc./box;
normcc = normcc/norm(eGf);
clear cc
toc;

box = box(winlen:end);
normcc = normcc(winlen:end);
toc;

%%
close all;
[pks,locs] = findpeaks(normcc,'MINPEAKDISTANCE',length(eGf),'MINPEAKHEIGHT',0.85);


%%
clear fc omega matches resnorm
thresh = 15;
ampRat = box(locs)./norm(eGf);
lI = (1:length(pks))'; %find(ampRat >= thresh);
disp([ampRat(lI) normcc(locs(lI))]);

%%
lfc = 1;
hfc = 45;
close all;
n = 1;
disp('looping');
matches = NaN(length(eGf),length(lI));
if diffFlag
    dorig = cumsum(detrend(zpkFilter(diff(dorig),1,-inf,100)));
end

for i = 1:length(lI)
    disp(ampRat(lI(i)))
    tfind = find(t>=t(locs(lI(i))),1);
    dfound = dorig(tfind:tfind+length(eGf)-1);
    matches(:,i) = detrend(dfound);
end
badI = find(sum(~isfinite(matches)));
matches(:,badI) = [];
maxAmpRMS = rms(matches)';
[~,minI] = min(maxAmpRMS);
eGf = matches(:,minI);
% eGf = normalizeWaveforms(nanmean(matches,2));
% eGf = eGf / norm(eGf);

for i = 1:length(lI)
    tfind = find(t>=t(locs(lI(i))),1);
    datestr(t(tfind))
    dfound = matches(:,i);
    
    Fs = 100;
    xfact = 1;
    eGf = resample(eGf,xfact,1);
    dfound = resample(dfound,xfact,1);
    
    %%
    if pmtmFlag
        [pegf,fegf] = pmtm(parzenwin(length(eGf)).*detrend(eGf),{2,'trace'},2^(nextpow2(length(eGf))+2),Fs*xfact);
        [pd,fd] = pmtm(parzenwin(length(dfound)).*detrend(dfound),{2,'trace'},2^(nextpow2(length(dfound))+2),Fs*xfact);
    else
        [pegf, fegf] = pwelch(taper(detrend(eGf),tw),parzenwin(64),32,2^(nextpow2(length(eGf))+2),Fs*xfact);
        [pd, fd] = pwelch(taper(detrend(dfound),tw),parzenwin(64),32,2^(nextpow2(length(dfound))+2),Fs*xfact);
    end
    fI2 = fegf >= lfc & fegf < hfc;
    
    %%
    mrs = pd./pegf;
    if curveFitFlag
        if boatwright
            gam = 2;
            slope = 2;
            ngam = slope*gam;
            f = @(a,xdata) a(1)./(1+(xdata./a(2)).^ngam).^(1/2);
        else
            gam = 1;
            slope = 2;
            ngam = slope*gam;
            f = @(a,xdata) a(1)./(1+(xdata./a(2)).^ngam).^(1/gam);
        end
        tic;
        [ahat,resnorm(i)] = lsqcurvefit(f,[max(mrs(fI2)) 8],fegf(fI2),mrs(fI2));
        toc;
        omega_ = ahat(1);
        fc_ = ahat(2);
        d = fegf(fI2);
        mrs = mrs(fI2);
        
    else
        tic;
        d = fegf(fI2);
        mrs = mrs(fI2);
        if boatwright
            gam = 2;
            slope = 2;
            ngam = slope*gam;
            fun = @(a)(a(1)./(1+(d./a(2)).^ngam).^(1/gam)) - mrs;
        else %brune
            gam = 1;
            slope = 2;
            ngam = slope*gam;
            fun = @(a)(a(1)./(1+(d./a(2)).^ngam).^(1/gam)) - mrs;
        end
        a0 = [max(mrs) 8];
        [ahat,resnorm(i)] = lsqnonlin(fun,a0);
        omega_ = ahat(1);
        fc_ = ahat(2);
        toc;
    end
    
    if fc_ > hfc || fc_ < lfc
        disp('skipping')
        continue
    else
        disp(n)
        disp([fc_ omega_])
        figure(1);
        hold on;
        plot(dfound/norm(dfound),'r')
        figure(2); hold on;
        loglog(d,pegf(fI2))
        hold on;
        loglog(fd(fI2),pd(fI2),'r')
        figure(3);
        hold on;
        brune = omega_./(1+(d./fc_).^ngam).^(1/2);
        loglog(d,mrs,'r');
        
        %plot best least squares fit (brune spectrum)
        loglog(d,brune,'k--')
        legend('data','best fit')
        xlabel('frequency')
        ylabel('relative $\Omega_0$')
        %title(['\Omega_0: ',num2str(round(omega)),'; f_c: ',num2str(round(fc*10)/10)])
        plot(fc_,omega_./(1+(fc_./fc_).^ngam).^(1/2),'p','markersize',15,'markeredgecolor','k','markerfacecolor','y')
        fc(n) = fc_;
        omega(n) = omega_;
        n = n+1;
    end
end

figure(1);
hold on;
plot(eGf/norm(eGf),'k','linewidth',2)
title({'black=eGf';['number of matches found: ',num2str(length(lI))]})
zoom on;

figure(2)
set(gca,'Xscale','log','Yscale','log')
xlabel('frequency')
ylabel('spectral amplitude')
title({'velocity spectra, instrument response not removed';'black=eGf'})
zoom on;

data = [fc' omega'];
data = sortrows(data);
fc = data(:,1);
omega = data(:,2);
%P = polyfit(log10(fc),log10(omega),1);
P = robustfit(log10(fc),log10(omega));
omega2 = P(2)*log10(fc)+P(1);
omega2 = 10.^omega2;

figure(3);
hold on;
zoom on;
loglog(fc,omega2,'b','linewidth',2)
title(['slope of blue line: ',num2str(P(2))]);
set(gca,'Xscale','log','Yscale','log')

for i = 1:length(lI)
    tmp = detrend(matches(:,i));
    matches(:,i) = tmp/norm(tmp);
end
figure(1);
hold on;
plot(nanmean(matches,2),'b','linewidth',2);
zoom on;

%%
% hold on;
% fI = fegf >= 8 & fegf < hfc;
% P = polyfit(log10(fegf(fI)),log10(mrs(fI)),1);
% disp(P)
% mrs_synth = P(1)*log10(fegf(fI)) + P(2);
% mrs_synth = 10.^mrs_synth;
% % loglog(fegf(fI),mrs_synth,'r','linewidth',2)
%
% [Cxy,F] = mscohere(eGf,dfound,[],[],[],100);
% figure; plot(F,Cxy)
% stf_ = ifft(fft(dfound)./fft(eGf));
% [pd, fd] = pmtm(taper(detrend(stf_),0.005),{4,'trace'},2^(nextpow2(length(stf_))+1),Fs*xfact);
% figure; loglog(fd,pd,'c')
% for i = 1:499
%     stf = getSTF(i/1000,Fs*xfact);
%     dfound2 = conv(stf,dfound);
%     ccoeff(i) = max(xcorr(dfound2(1:length(eGf)),eGf,'coeff'));
% end