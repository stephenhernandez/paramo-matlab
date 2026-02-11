clear; close all; clc;
cd ~/research/now/sangay/

load saga_drumbeat_data.mat
[pxxUnfiltered,f] = pmtm(dUnfiltered,2.5,2^15,100,'unity');
pxxSmall = pxxUnfiltered(:,amps < 5000);
close all; 
figure(); ax = gca; loglog(f,pxxUnfiltered(:,amps>=5000),'linewidth',0.1,'Color',ax.ColorOrder(1,:)); zoom on; hold on; 
loglog(f,pxxUnfiltered(:,amps<5000),'linewidth',0.1,'Color',ax.ColorOrder(2,:)); zoom on;

%
averageEGF_spectra = nanmedian(pxxSmall,2);
hold on; loglog(f,averageEGF_spectra,'k','linewidth',4); zoom on;
spectralRatios = pxxUnfiltered./averageEGF_spectra;
figure(); ax = gca; loglog(f,spectralRatios(:,amps>=5000),'linewidth',0.1,'Color',ax.ColorOrder(1,:)); zoom on; hold on; loglog(f,spectralRatios(:,amps<5000),'linewidth',0.1,'Color',ax.ColorOrder(2,:)); zoom on;
hold on; loglog(f,nanmean(spectralRatios(:,amps>=5000),2),'k','linewidth',2); zoom on; hold on; loglog(f,nanmean(spectralRatios(:,amps<5000),2),'k','linewidth',2); zoom on;
fI = f >= 0.7 & f <= 30;

%%
figure();
semilogy(t2(amps >= 5000),max(spectralRatios(:,amps>=5000)),'.'); zoom on;

%
aI = amps >= 4000;
aI = find(aI);
la = length(aI); d = double(f(fI)); spectralRatios = double(spectralRatios); fc = NaN(la,1); omega = fc; slopey = fc; resnorm = fc; for i = 1:la
    tic;
    mrs = spectralRatios(fI,aI(i));
    %
    gam = 2;
    slope = 2;
    ngam = slope*gam;disp('i am here');
    fun = @(a,xdata) a(1)./(1+(xdata./a(2)).^ngam).^(1/gam);
    [ahat,resnorm(i)] = lsqcurvefit(fun,[max(mrs) 10],d,mrs);
    omega_ = ahat(1); omega(i) = omega_;
    fc_ = ahat(2); fc(i) = fc_;
    toc;
end

%%
figure(); 
i2 = fc > 1/20 & omega >= 10; loglog(1./fc(i2),omega(i2),'.'); zoom on;
b_ = robustfit(log10(1./fc(i2)),log10(omega(i2)));
b_ = flipud(b_);
yq = polyval(b_,log10(1./fc(i2)));
figure(4); hold on; ll = loglog(1./fc(i2),10.^yq,'linewidth',2); ll.Color(4) = 0.6; title(['slope = ',num2str(b_(1))]);
grid on;
