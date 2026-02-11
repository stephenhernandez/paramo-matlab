%plot_best_solution
%function plot_best_solution(ref,t0_,th_)
close all
ecuadorFlag = true;
% t0_ = 6;
% th_ = 1;
% i_ = TH <= 40;
% error_tmp = error_tmp(:,i_);
% T0 = T0(i_);
% TH = TH(i_);
% dc = dc(:,i_);

% if dcFlag
%     %error = error_tmp.*dc/100;
%     error = max(error_tmp.*dc)/100;
% else
%     error = max(error_tmp);
% end

maxerror = max(error);
if ~ecuadorFlag
    mI = find(error_tmp == maxerror);
else
    mI = find(error == maxerror);
end
%mI = find(T0 == t0_ & TH == th_);

mI = mI(end);
t0_ = T0(mI);
th_ = TH(mI);

%%
tI = find(T0 == t0_);
tI2 = find(TH == th_);
%load ~/igdata/Top16.mat;
%load data/mat.mat;
close all;
if threeDFlag
    depthI = depth == depth(maxErrI);
end

if dcFlag
    figure('units','normalized','outerposition',[0 0 1 1]);
    if threeDFlag
        scatter(lons(depthI),lats(depthI),80,error_tmp(depthI,mI).*dc(depthI,mI)/100,'filled'); h = colorbar;  zoom on; grid on;%caxis([0 100]);
        [~,locI] = max(error_tmp(:,mI).*dc(:,mI));
    else
        scatter(lons,lats,80,error_tmp(:,mI).*dc(:,mI)/100,'filled'); h = colorbar;  zoom on; grid on;%caxis([0 100]);
        [~,locI] = max(error_tmp(:,mI).*dc(:,mI));
    end
else
    figure('units','normalized','outerposition',[0 0 1 1]);
    if threeDFlag
        scatter(lons(depthI),lats(depthI),80,error_tmp(depthI,mI),'filled'); h = colorbar;  zoom on; grid on;%caxis([0 100]);
        [~,locI] = max(error_tmp(:,mI));
    else
        scatter(lons,lats,80,error_tmp(:,mI),'filled'); h = colorbar;  zoom on; grid on;%caxis([0 100]);
        [~,locI] = max(error_tmp(:,mI));
    end
end

ylabel(h,'Variance Reduction');
hold on;
geoshow('~/igdata/fallas/fallas2016.shp','Color',[0.5 0.5 0.5]);
hold on; 
geoshow('~/igdata/ZonaUrbana/ZonaUrbana.shp');
xlabel('Longitude');
ylabel('Latitude');
nCmpsUsed = sum(keepI);
plot(lonEC,latEC,'k','linewidth',2)

n0 = find(nCmpsUsed == 0);
n1 = find(nCmpsUsed == 1);
n2 = find(nCmpsUsed == 2);
n3 = find(nCmpsUsed == 3);

if ~isempty(n0)
    plot(stlo(n0),stla(n0),'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k','markersize',10);  zoom on; grid on;
    for kk = 1:length(n0)
        text(stlo(n0(kk)),stla(n0(kk)),stnms{n0(kk)},'FontSize',10)
    end
end
if ~isempty(n1)
    plot(stlo(n1),stla(n1),'s','markerfacecolor','w','markeredgecolor','k','markersize',10);
    for kk = 1:length(n1)
        text(stlo(n1(kk)),stla(n1(kk)),stnms{n1(kk)},'FontSize',10)
    end
end
if ~isempty(n2)
    plot(stlo(n2),stla(n2),'d','markerfacecolor','w','markeredgecolor','k','markersize',10);
    for kk = 1:length(n2)
        text(stlo(n2(kk)),stla(n2(kk)),stnms{n2(kk)},'FontSize',10)
    end
end
if ~isempty(n3)
    plot(stlo(n3),stla(n3),'v','markerfacecolor','w','markeredgecolor','k','markersize',10);
    for kk = 1:length(n3)
        text(stlo(n3(kk)),stla(n3(kk)),stnms{n3(kk)},'FontSize',10)
    end
end

%hold on; plot(mat(:,1),mat(:,2),'k--','linewidth',2);

axis equal;
axis([minlon maxlon minlat maxlat]);
ldur = length(dur);
lBig = ldur*durations(1,mI);
if flag3
    lBig = lBig*3;
end
plot_waveform(synthBig(1:lBig,mI),obsBig(1:lBig,mI),durations(1:ldur,mI),dists(1:ldur,mI),stnms(1:ldur),200,2,fs,flag3,keepI);
%plotWaveform(synthSV(1:lBig,mI),obsBig(1:lBig,mI),durations(1:ldur,mI),dists(1:ldur,mI),stnms(1:ldur),50,1.5,fs,flag3)

disp(['Best Lon.: ',num2str(lons(locI))]);
disp(['Best Lat.: ',num2str(lats(locI))]);
disp(['Best Depth.: ',num2str(depth(locI))]);
disp(['Best Origin Time: ',datestr(ref+(T0(mI)/86400)),', (',num2str(T0(mI)),')']);
disp(['Best Half-Duration: ',num2str(TH(mI))]);
disp(['Best Mw: ',num2str(mw(mI))]);
disp(['Best Variance Reduction: ',num2str(error_tmp(locI,mI))]);
disp(['Best DC: ',num2str(dc(locI,mI))]);
disp(['Best DC-Scaled Variance Reduction: ',num2str(error_tmp(locI,mI).*dc(locI,mI)/100)]);

figure(1);
grid on;
if ~ecuadorFlag
    title(['Lat: ' num2str(lats(locI)), ', Lon: ' num2str(lons(locI)),...
        ', Depth: ' num2str(depth(locI)),', $M_W$: ', num2str(mw(mI)),', wVR: ',num2str(round(error_tmp(mI)))]);
else
    title(['Lat: ' num2str(lats(locI)), ', Lon: ' num2str(lons(locI)),...
        ', Depth: ' num2str(depth(locI)),', $M_W$: ', num2str(mw(mI)),', wVR: ',num2str(round(error(mI)))]);
end

G = Gbig_tmp;
gtg = G'*G;
[V,~] = eig(gtg);

if threeDFlag
    lI = lons == lons(locI) & lats == lats(locI);
    figure('units','normalized','outerposition',[0 0 1 1]);
    if dcFlag
        plot(error_tmp(lI,mI).*dc(lI,mI)/100,-depth(lI),'k.-'); zoom on; grid on;
    else
        plot(error_tmp(lI,mI),-depth(lI),'k.-'); zoom on; grid on;
    end
    xlabel('Variance Reduction [$\%$]');
    ylabel('Depth [km.]')
end

%% plot concatenated solution
figure('units','normalized','outerposition',[0 0 1 1]);
plot(obsBig(1:lBig,mI),'k','linewidth',2); zoom on; grid on;
hold on;
plot(synthBig(1:lBig,mI),'linewidth',2); zoom on; grid on;
legend('Obs.','Synth');

%%
figure('units','normalized','outerposition',[0 0 1 1]);
if ecuadorFlag
    scatter(TH(tI),error(tI),200,mw(tI),'filled'); h = colorbar; zoom on; grid on;
else
    scatter(TH(tI),error_tmp(tI),200,mw(tI),'filled'); h = colorbar; zoom on; grid on;
end
ylabel(h,'Moment Magnitude');
xlabel(['Half Duration (Origin Time fixed at ',datestr(ref+(t0_/86400)),')'])
ylabel('Variance Reduction (weighted by DC)');

figure('units','normalized','outerposition',[0 0 1 1]);
if ~ecuadorFlag
    scatter(origt+seconds(T0(tI2)),error_tmp(tI2),200,mw(tI2),'filled'); h = colorbar; zoom on; grid on;
else
    scatter(origt+seconds(T0(tI2)),error(tI2),200,mw(tI2),'filled'); h = colorbar; zoom on; grid on;
end
ylabel(h,'Moment Magnitude');
xlabel(['Origin Time (Half Duration fixed at ',num2str(th_),')'])
ylabel('Variance Reduction (weighted by DC)');
datetick('x','keeplimits');

figure('units','normalized','outerposition',[0 0 1 1]);
if ecuadorFlag
    scatter(origt+seconds(T0),TH,80,error,'filled'); h = colorbar; zoom on; grid on;
else
    scatter(origt+seconds(T0),TH,80,error_tmp,'filled'); h = colorbar; zoom on; grid on;
end

ylabel(h,'variance reduction');
xlabel('Origin Time')
ylabel('Half-duration');
datetick('x','keeplimits');
hold on; plot(origt+seconds(T0(mI)),TH(mI),'p','markersize',20,'markeredgecolor','w','markerfacecolor','k');  zoom on; grid on;
axis tight


% % if threeDFlag
% %     lI = lons == lons(locI) & lats == lats(locI);
% %     figure('units','normalized','outerposition',[0 0 1 1]);
% %     if dcFlag
% %         plot(error_tmp(lI,mI).*dc(lI,mI)/100,-depth(lI),'k.-')
% %     else
% %         plot(error_tmp(lI,mI),-depth(lI),'k.-')
% %     end
% %     xlabel('Variance Reduction [$\%$]');
% %     ylabel('Depth [km.]')
% % end

[strike,dip,rake] = mt2sdr(mBig(:,mI)');
[strike_a,dip_a,rake_a] = auxplane([strike, dip, rake]);
% [strike,dip,rake] = mt2sdr(mBig(:,mI)');
% [strike,dip,rake] = auxplane([strike,90-dip, rake]);
% [strike_a,dip_a,rake_a] = auxplane([strike,dip, rake]);
disp('nodal plane 1')
disp(round([strike dip rake]))
disp('nodal plane 2')
disp(round([strike_a dip_a rake_a]))
