function [mc,b,bmin,bmax,alphaHat,alphaMOM,betaHat,betaMOM,mtest,momCounts,mom1] = mcGamma(mag,nboot,plotFlag,dM)
%define default outputs
if nargin < 4; dM = 0.01; end
if nargin < 3; plotFlag = false; end
disp(' ');

% clear; close all; clc;
% cd ~/research/now/pedernales_forecast/
% load final_aftershock_catalog.mat

% dM = 0.01;
% nboot = 2e2;
% t2 = t(goodI);
% tRef = t2(newmeanmag2 == max(newmeanmag2));
% tI = (t2 > tRef);
% newmeanmag2 = round(newmeanmag/dM)*dM;
%
% mag = newmeanmag2(goodI);
% mag = mag(tI);
% mag = mag(1:100);
% plotFlag = true;


%%
dm2 = dM/2;
mag = round(mag/dM)*dM;
pr2 = sort(mag);
N = length(mag);
C = log10(exp(1));
%figure(); N = length(mag); plot(sort(mag),(N - (1:N)')/N,'o'); zoom on; grid on;

%%
minMag = min(mag)+dM;
maxMag = max(mag)+dM;
%maxmtest = max([minMag+0.8 maxMag-0.8]);
%maxmtest = max([minMag+1.5 maxMag-1.5]);
maxmtest = max([minMag+1 maxMag-1]);

%%
mtest = unique(mag); %round(mag*10)/10);
mtI = mtest <= maxmtest;
mtest = sort(mtest(mtI));

Ntest = NaN(sum(mtI),1);
for i = 1:length(mtest)
    Ntest(i) = sum(mag >= mtest(i));
end
Nmin = 50;
mtI = Ntest >= Nmin;
mtest = sort(mtest(mtI));

% Nmin = 50;
% if N < Nmin
%     disp('Aborting');
%     return
% end

Nmtest = length(mtest);
mom1 = NaN(Nmtest,nboot);
mom2 = mom1;
mom3 = mom1;

%% preallocate
alphaHat = mom1;
betaHat = mom1;
betaMOM = betaHat;
alphaMOM = betaHat;
gammaMOM = mom1;

%%
N
bootI = randi(N,N,nboot);
mscramble = mag(bootI);
momCounts = zeros(1,Nmtest);
for j = 1:nboot
    msub = mscramble(:,j);
    for i = 1:Nmtest
        mI = msub >= mtest(i);
        minX = min(msub);
        X = msub(mI) - minX + dm2;
        %newMinX = min(X);

        % moments
        xbar = mean(X);
        mom1_ = xbar;                   % first raw moment
        mom2_ = mean((X - xbar).^2);    % second central moment
        mom3_ = mean((X - xbar).^3);    % third central moment
        
        %disp([mom1_ mom2_ mom3_]);
        mom1(i,j) = mom1_;
        mom2(i,j) = mom2_;
        mom3(i,j) = mom3_;
        
        % parameter estimation using method of moments
        alphaMOM(i,j) = 4*(mom2_^3)/(mom3_^2);
        alphaHat(i,j) = alphaMOM(i,j);
        betaMOM(i,j) = 2*mom2_/mom3_;
        betaHat(i,j) = 2*mom2_/mom3_;
        gamma_ = mom1_ - (2*(mom2_^2)/mom3_);
        gammaMOM(i,j) = gamma_;
        
        %%
        Y = X - gamma_;
        mI = Y > 0;
        Y = Y(mI);

        logY = log(Y);
        meanlogY = mean(logY);
        Ybar = mean(Y);
        logYbar = log(Ybar);
        minY = min(Y);

        MM = logYbar - meanlogY;
        k_ = (3 - MM + sqrt(((MM-3)^2)+24*MM))/(12*MM);
        %disp([k_ alphaMOM(i,j)])
        f = @(a) log(a) - psi(a) - MM;
        
        % here i try to estimate mle
        try
            if isfinite(k_)
                alphaHat_ = fzero(f,k_);%alphaMOM(i,j));
            else
                alphaHat_ = fzero(f,1); %alphaMOM(i,j));
            end
            alphaHat(i,j) = alphaHat_;
            %betaHat(i,j) = alphaHat(i,j)./mom1(i,j);
            %betaHat(i,j) = 1./Ybar; %./(alphaHat(i,j));
            betaHat(i,j) = exp(psi(alphaHat_) - meanlogY);
            %disp([alphaHat_ alphaMOM(i,j) betaHat(i,j) betaMOM(i,j) gamma_ gammaMOM(i,j) minY newMinX]);
        catch ME
            disp('didnt work, skipping');
            [k_ alphaMOM(i,j)]
            fprintf('%f %f %f %f %f %f %f %f\n',alphaMOM(i,j),betaMOM(i,j),gamma_,gammaMOM(i,j),minY,meanlogY,logYbar,Ybar);
            rethrow(ME)
            fprintf('%f %f %f %f\n',alphaMOM(i,j),betaMOM(i,j),betaHat(i,j),alphaHat_);
            alphaHat(i,j) = NaN;
            
%             try
%                 alphaHat(i,j) = fzero(fun,0.5);
%                 %betaHat(i,j) = alphaHat(i,j)./mom1(i,j);
%                 %betaHat(i,j) = 1./Ybar; %./(alphaHat(i,j));
%                 betaHat(i,j) = 1./exp(meanlogY - psi(alphaHat(i,j)));
%             catch
%                 try
%                     alphaHat(i,j) = fzero(fun,1);
%                     %betaHat(i,j) = alphaHat(i,j)./mom1(i,j);
%                     %betaHat(i,j) = 1./Ybar; %./(alphaHat(i,j));
%                     betaHat(i,j) = 1./exp(meanlogY - psi(alphaHat(i,j)));
%                 catch
%                     try
%                         alphaHat(i,j) = fzero(fun,2);
%                         %betaHat(i,j) = alphaHat(i,j)./mom1(i,j);
%                         %betaHat(i,j) = 1./Ybar; %./(alphaHat(i,j));
%                         betaHat(i,j) = 1./exp(meanlogY - psi(alphaHat(i,j)));
%                     catch
%                         try
%                             alphaHat(i,j) = fzero(fun,3);
%                             %betaHat(i,j) = alphaHat(i,j)./mom1(i,j);
%                             %betaHat(i,j) = 1./Ybar; %./(alphaHat(i,j));
%                             betaHat(i,j) = 1./exp(meanlogY - psi(alphaHat(i,j)));
%                         catch
%                             %disp('mle not possible. using method of moments.')
%                             %disp([mtest(i) j])
%                             alphaHat(i,j) = alphaHat(i,j);
%                             %betaHat(i,j) = betaHat(i,j);
%                             %betaHat(i,j) = 1./Ybar; %./(alphaHat(i,j));
%                             betaHat(i,j) = 1./exp(meanlogY - psi(alphaHat(i,j)));
%                             momCounts(i) = momCounts(i) + 1;
%                         end
%                     end
%                 end
%             end
        end
    end
end

%%
centroid_alpha = median(alphaHat,2,"omitnan");
centroid_beta = median(betaHat,2,"omitnan");
centroid_gamma = median(gammaMOM,2,"omitnan");
centroid_alphaGamma = median(alphaHat - gammaMOM - 1,2,"omitnan");

ciAlpha = 1;
alphaLB = prctile(alphaHat',ciAlpha/2);
alphaUB = prctile(alphaHat',100 - (ciAlpha/2));

figure(); plot(alphaHat(1,:),'.'); zoom on; grid on;
figure(); plot(betaHat(1,:),'.'); zoom on; grid on;
betaLB = prctile(betaHat',ciAlpha/2);
betaUB = prctile(betaHat',100-ciAlpha/2);
gammaLB = prctile(gammaMOM',ciAlpha/2);
gammaUB = prctile(gammaMOM',100-ciAlpha/2);

alphaGammaLB = prctile((alphaHat - gammaMOM - 1)',ciAlpha/2);
alphaGammaUB = prctile((alphaHat - gammaMOM - 1)',100-ciAlpha/2);

%%
chooseMetric = abs(centroid_alpha - 1).*abs(centroid_alphaGamma);
mcI = find(alphaLB <= 1 & alphaUB >= 1,1);% | (alphaGammaLB <= 0 & alphaGammaUB >= 0));
if isempty(mcI)
    mcI = find(alphaGammaLB <= 0 & alphaGammaUB >= 0,1);
end

if isempty(mcI)
    disp('no proper mc found')
    mcI = find(chooseMetric == min(chooseMetric));
elseif length(mcI) > 1
    disp(['candidate mcs: ',num2str(mtest(mcI)')]);
    chooseMetric_ = abs(chooseMetric(mcI));
    [~,min_] = min(chooseMetric_);
    mcI = mcI(min_); %closest to zero
end

%%
mc = mtest(mcI);
b = C.*centroid_beta(mcI);
b = round(b*100)/100;
bdist = C.*betaHat(mcI,:);
y = prctile(bdist,[ciAlpha/2 100-ciAlpha/2]);
bmin = y(1);
bmax = y(2);
bmin = round(bmin*100)/100;
bmax = round(bmax*100)/100;

%% plot if desired
if plotFlag
    pentMarkerSize = 25;
    figure('units','normalized','outerposition',[0 0 1 0.75]);
    ax(1) = subplot(3,2,[2 4 6]);
    curveGR = (N - (0:N-1)');
    semilogy(pr2,curveGR,'o-'); grid on; zoom on; hold on;
    mcpr = find(pr2 >= mtest(mcI),1);
    pp = plot(pr2(mcpr)*[1 1],[0.1 10^(ceil(log10(max(curveGR))))],'k-','linewidth',5); pp.Color(4) = 0.25;
    pp = plot([min(pr2) max(pr2)],curveGR(mcpr)*[1 1],'k-','linewidth',5); pp.Color(4) = 0.25;
    plot(pr2(mcpr),curveGR(mcpr),'p','MarkerSize',pentMarkerSize,'linewidth',5);
    xsynth = (pr2(mcpr):dM:maxMag)';
    ysynth = -b.*(xsynth - pr2(mcpr))+log10(curveGR(mcpr));
    ysynth1 = -bmax.*(xsynth - pr2(mcpr))+log10(curveGR(mcpr));
    ysynth2 = -bmin.*(xsynth - pr2(mcpr))+log10(curveGR(mcpr));
    
    semilogy(xsynth,10.^ysynth,'.');
    semilogy(xsynth,10.^ysynth1,'-','color',[0.5 0.5 0.5]);
    semilogy(xsynth,10.^ysynth2,'-','color',[0.5 0.5 0.5]);
    
    title(['$M_{C}$: ',num2str(mtest(mcI)),', $b \textendash value$: ',num2str(b),', $b_{min}$: ',num2str(bmin),', $b_{max}$: ',num2str(bmax)]);
    
    ax(2) = subplot(3,2,1);
    hold on;
    pp = plot([minMag-2*dM maxMag],[1 1],'k-','linewidth',5); pp.Color(4) = 0.25;
    plot(mtest,alphaLB,'.-','color',[0.5 0.5 0.5],'linewidth',2);
    plot(mtest,alphaUB,'.-','color',[0.5 0.5 0.5],'linewidth',2);
    plot(mtest,centroid_alpha,'.-','linewidth',2); grid on; zoom on;
    plot(mtest(mcI),centroid_alpha(mcI),'p','MarkerSize',pentMarkerSize,'linewidth',5);
    ylim([0 4]);
    ylabel('$\hat{\alpha}$');
    
    ax(3) = subplot(3,2,3);
    hold on;
    pp = plot([minMag-2*dM maxMag],0*[1 1],'k-','linewidth',5); pp.Color(4) = 0.25;
    plot(mtest,gammaLB,'.-','color',[0.5 0.5 0.5],'linewidth',2);
    plot(mtest,gammaUB,'.-','color',[0.5 0.5 0.5],'linewidth',2);
    plot(mtest,centroid_gamma,'.-','linewidth',2); grid on; zoom on;
    plot(mtest(mcI),centroid_gamma(mcI),'p','MarkerSize',pentMarkerSize,'linewidth',5);
    ylim([-3 3]);
    ylabel('$\hat{\gamma}$');
    linkaxes(ax,'x');
    zoom on;
    
    ax(4) = subplot(3,2,5);
    hold on;
    pp = plot([minMag-2*dM maxMag],0*[1 1],'k-','linewidth',5); pp.Color(4) = 0.25;
    plot(mtest,alphaGammaLB,'.-','color',[0.5 0.5 0.5],'linewidth',2);
    plot(mtest,alphaGammaUB,'.-','color',[0.5 0.5 0.5],'linewidth',2);
    plot(mtest,centroid_alphaGamma,'.-','linewidth',2); grid on; zoom on;
    plot(mtest(mcI),centroid_alphaGamma(mcI),'p','MarkerSize',pentMarkerSize,'linewidth',5);
    ylim([-1 3]);
    ylabel('$\hat{\alpha} - \hat{\gamma} - 1$');
    linkaxes(ax,'x');
    zoom on;
    
    figure('units','normalized','outerposition',[0.2 0.2 0.75 0.75]);
    plot(mtest,C.*betaHat,'-','color',[0.5 0.5 0.5],'linewidth',0.1);
    hold on;
    plot(mtest,C.*centroid_beta,'s-','linewidth',2);
    plot(mtest,C.*betaLB,'.-','color','k','linewidth',2);
    plot(mtest,C.*betaUB,'.-','color','k','linewidth',2);
    plot(mtest(mcI),b,'p','MarkerSize',pentMarkerSize,'linewidth',5);
    title(['$M_{C}$: ',num2str(mtest(mcI)),', $b \textendash value$: ',num2str(b),', $b_{min}$: ',num2str(bmin),', $b_{max}$: ',num2str(bmax)]);
    ylim([0 3]);
    ylabel('b-value'); % (median of several bootstraps)');
    zoom on;
    grid on;
end
