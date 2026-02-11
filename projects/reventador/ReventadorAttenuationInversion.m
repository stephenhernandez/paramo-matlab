%ReventadorAttenuationInversion
clear; close all;
cd ~/research/now/reventador/
load correctedAmplitudes.mat

%%
boniCorrection = 1201/749;
col = 18; %<-- 16 uses undocumented spectral ratio correction
Rcorr_z2p = [allCASCAmps(:,col) boniCorrection*allBONIAmps(:,col) allANTSAmps(:,col) allANTGAmps(:,col)];
goodEventsI = sum(isfinite(Rcorr_z2p)&Rcorr_z2p>=10,2);
fourI = goodEventsI == 4;
fourI2 = find(fourI);

subI = sort(fourI2(randi(length(fourI2),1.2e4,1)));

%
d_ = d_*1e-3;
% figure(); hold on; 
% for i = 1:4
%     loglog(d_(i),Rcorr_z2p(subI,i),'.'); zoom on; grid on;
% end
% ax = gca;
% ax.YScale = 'log';
% ax.XScale = 'log';

G_ = full(Gvdcc(4));
G_ = G_(1:end-1,:);
G_ = [getDD(d_) getDD(log10(d_)) G_];
%G_ = [getDD(log10(d_)) G_];

mbest_ = [];
rdum = logspace(0,3,1001);
tic;
for i = 1:21
    col = i;
    disp([i; median(allCASCAmps(fourI,col)./(boniCorrection*allBONIAmps(fourI,col)),"omitnan");...
        median(allCASCAmps(fourI,col)./(allANTSAmps(fourI,col)),"omitnan");...
        median(allCASCAmps(fourI,col)./(allANTGAmps(fourI,col)),"omitnan");...
        median((boniCorrection*allBONIAmps(fourI,col))./allANTSAmps(fourI,col),"omitnan");...
        median((boniCorrection*allBONIAmps(fourI,col))./allANTGAmps(fourI,col),"omitnan");...
        median((allANTSAmps(fourI,col))./allANTGAmps(fourI,col),"omitnan")]);

    ampDum = [allCASCAmps(:,col) boniCorrection*allBONIAmps(:,col) allANTSAmps(:,col) allANTGAmps(:,col)];
    G = [];
    d = [];
    for j = 1:length(subI)
        amps_ = [allCASCAmps(subI(j),col) boniCorrection*allBONIAmps(subI(j),col) allANTSAmps(subI(j),col) allANTGAmps(subI(j),col)];
        G = [G; G_];
        d = [d; getDD(log10(amps_'))];
        %disp(j);
    end

    badD = find(~isfinite(d));
    G(badD,:) = [];
    d(badD) = [];
    G = [G; 0 0 1 1 1 1];
    d = [d; 0];

    %
    %mbest = ((G'*G)^(-1))*G'*d;
    mbest = pinv(G)*d;
    mbest_ = [mbest_ mbest];
    disp(mbest_);
    figure(1);
    hold on;
    att = mbest(1).*rdum+mbest(2)*log10(rdum);
    if any(~isfinite(mbest(1:2)))
        continue;
    end

    plotSymbol = '.';
    if i < 15
        plotSymbol = 'p';
    elseif i < 22
        plotSymbol = 'p';
    end
    semilogx(rdum,-att,plotSymbol); zoom on; grid on;

    Mcorr = [allCASCAmps(fourI,col) boniCorrection*allBONIAmps(fourI,col) allANTSAmps(fourI,col) allANTGAmps(fourI,col)];
    for j = 1:4
        Mcorr(:,j) = log10(Mcorr(:,j))-(mbest(1)*d_(j)+mbest(2)*log10(d_(j))+mbest(j+2));
        %Mcorr(:,j) = log10(Rcorr_z2p(:,j))-(mbest(1)*log10(d_(j))+mbest(j+1));
    end
    fprintf('%d, mean mad: %g\n',i,mean(mad(Mcorr,1,2)));
    clear Mcorr;
    tic;
end
toc;

%%
figure('units','normalized','outerposition',[0 0 1 1/2]); hold on;
for i = 1:3
    ref = Rcorr_z2p(fourI,i);
    for j = i+1:4
        H = histogram(ref./Rcorr_z2p(fourI,j),logspace(-1,2,1501)); zoom on; grid on; ax = gca; ax.XScale = 'log';
        H.EdgeColor = 'none';
        H.FaceAlpha = 0.5;
    end
end

%%
G = []; 
d = [];
for j = 1:length(subI)
    amps_ = Rcorr_z2p(subI(j),:);
    G = [G; G_];
    d = [d; getDD(log10(amps_'))];
    disp(j);
end

G = [G; 0 0 1 1 1 1];
d = [d; 0];
rdum = logspace(0,3,1001);

%
%mbest = ((G'*G)^(-1))*G'*d;
mbest = pinv(G)*d;
fprintf('mean correction: %f\n',mean(mbest(3:end)));

%%
for i = 1:21
    Rcorr_z2p = [allCASCAmps(fourI,i) boniCorrection*allBONIAmps(fourI,i) allANTSAmps(fourI,i) allANTGAmps(fourI,i)];
    Mcorr = Rcorr_z2p;
    for j = 1:4
        Mcorr(:,j) = log10(Rcorr_z2p(:,j))-(mbest(1)*d_(j)+mbest(2)*log10(d_(j))+mbest(j+2));
        %Mcorr(:,j) = log10(Rcorr_z2p(:,j))-(mbest(1)*log10(d_(j))+mbest(j+1));
    end
    fprintf('%d, mean mad: %g\n',i,mean(mad(Mcorr,1,2)));
end

%%
i = 18; %16 is the best but uses undocumented spectral ratio correction, 20 aint bad
Rcorr_z2p = [allCASCAmps(:,i) boniCorrection*allBONIAmps(:,i) allANTSAmps(:,i) allANTGAmps(:,i)];
Mcorr = Rcorr_z2p;
for j = 1:4
    Mcorr(:,j) = log10(Rcorr_z2p(:,j))-(mbest(1)*d_(j)+mbest(2)*log10(d_(j))+mbest(j+2));
end
M = median(Mcorr,2,"omitnan");

%%
nPerEvent = sum(goodEventsI,2);

writeFlag = false;
if writeFlag
    t = datenum(t_);
    t = t-693960;
    formatSpec = '%f %f %f %d';
    str = compose(formatSpec,t,median(Mcorr,2,"omitnan"),mean(Mcorr,2,"omitnan"),nPerEvent);
    str = string(str);
    fileID = fopen('~/ReveMagnitudes_v1.txt','w');
    fprintf(fileID,'%s\n',str);
    fclose(fileID);
end

