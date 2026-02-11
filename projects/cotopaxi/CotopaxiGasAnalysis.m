clear; close all; clc;
%cd ~/research/now/cotopaxi/;
cd ~/igdata;
nDays = 30;
T = readtable('CotopaxiData.xlsx');
R = T.Var2;
RS = T.Var4;
C = T.Var6;
S = T.Var8;
G = [R RS C S];

G(G==0) = NaN;
N = sum(isfinite(G),2);
maxgas = medfiltSH(nanmax(G,[],2,'omitnan'),nDays);
medgas = medfiltSH(nanmedian(G,2),nDays);
W = T.Var14;
t = T.Var1;

%
figure(); 
SS = scatter(maxgas,medfiltSH(W,nDays),100,datenum(t),'filled'); 
zoom on; grid on; Cbar = colorbar; SS.MarkerFaceAlpha = 0.5; 
Wax = gca; Wax.XScale = 'log';
xlabel('median max gas [7-day window]');
ylabel('median windspeed [7-day window]');

%
figure('units','normalized','outerposition',[1/2 0 1/4 1]);
ax3(1) = subplot(311); semilogy(t,maxgas,'.'); 
zoom on; grid on;
title('7-day median of daily Max Gas value');
ax3(2) = subplot(312); semilogy(t,medgas,'.');
zoom on; grid on; title('7-day median of daily median gas value');
ax3(3) = subplot(313); plot(t,medfiltSH(W,nDays),'.'); zoom on; grid on; 
title('7-day median of daily average wind speed'); 
linkaxes(ax3,'x');

%
figure('units','normalized','outerposition',[0 0 1/4 1]); 
ys = (2015:2021)';
lys = length(ys);
for i = 1:lys
    ax(i) = subplot(lys,1,i); 
    tI = t>=datetime(ys(i),01,01) & t< datetime(ys(i)+1,01,01); 
    ll = plot(ax(i),t(tI),medgas(tI),'linewidth',4); 
    ll.Color(4) = 0.5; zoom on; grid on;
    ax(i).YScale = 'log'; 
end
linkaxes(ax,'y');
title(ax(1),'median gas: blue; max gas: red');

figure(3); hold on; 
ys = (2015:2021)'; 
lys = length(ys); 
for i = 1:lys
    ax(i) = subplot(lys,1,i); hold on; tI = t>=datetime(ys(i),01,01) & t< datetime(ys(i)+1,01,01); 
    ll = plot(ax(i),t(tI),maxgas(tI),'linewidth',4); 
    ll.Color(4) = 0.5; zoom on; grid on; ax(i).YScale = 'log'; 
end
linkaxes(ax,'y');

%figure 4
figure('units','normalized','outerposition',[3/4 0 1/4 1]); 
ys = (2015:2021)'; 
lys = length(ys); 
for i = 1:lys
    ax(i) = subplot(lys,1,i); 
    tI = t>=datetime(ys(i),01,01) & t< datetime(ys(i)+1,01,01); 
    ll = plot(ax(i),t(tI),W(tI),'.'); 
    zoom on; grid on; ax(i).YScale = 'lin'; 
    ylabel(ax(i),'$m/s$')
end
linkaxes(ax,'y');
title(ax(1),'Windpeeds');

% figure 5
figure();
axSmooth(1) = subplot(211);
plot(t,maxgas./medgas,'.'); 
zoom on; grid on; title('Smoothed Max / Smoothed Median');

figure(); plot(t,N,'.'); zoom on;

SensorNames = ["Refugio";"Refugio Sur";"Cami";"San Joaquin"];

Gs = G; 
for i = 1:4
    G_ = G(:,i); 
    G_ = medfiltSH(G_,nDays);
    Gs(:,i) = G_; 
end

figure(); 
hold on; n = 0; 
for i = 1:3 
    n = n+1; 
    subplot(3,2,n); t1 = t <= datetime(2017,04,01); 
    ratioGS = Gs(t1,4)./Gs(t1,i); 
    
    histogram(ratioGS,logspace(-1,2,101)); 
    ax_ = gca; ax_.XScale = 'log'; 
    zoom on; title(['Mean Ratio: ',num2str(nanmedian(ratioGS))]); 
    legend(strcat(SensorNames(4),'/',SensorNames(i))); 

    n = n+1; subplot(3,2,n); 
    t1 = t >= datetime(2018,10,01); 
    ratioGS = Gs(t1,4)./Gs(t1,i); 
    histogram(ratioGS,logspace(-1,2,101)); zoom on; 
    title(['Mean Ratio: ',num2str(nanmedian(ratioGS))]); 
    legend(strcat(SensorNames(4),'/',SensorNames(i))); ax_ = gca; ax_.XScale = 'log'; 
end
%ax4 = gca; ax4.YScale = 'log';
sgtitle('Left Column: t $\leq$ April 2017; Right Column: t $\geq$ October 2018','FontSize',24,'Interpreter','latex');

figure(); 
hold on; 
for i = 1:4
    G_ = G(:,i); 
    ll = plot(t,medfiltSH(G_,nDays),'linewidth',4); 
    ll.Color(4) = 0.75; zoom on; grid on; 
end
legend(SensorNames); sax = gca; sax.YScale = 'log';
title('7-day moving median');

figure(); 
slats = [-0.66373; -0.73462; -0.67980; -0.69200]; 
slons = [-78.44077; -78.47780; -78.50642; -78.58108]; 
plot(slons,slats,'v'); axis equal; 
text(slons,slats,SensorNames,'FontSize',20); grid on;

%% "corrected" version
t1 = t <= datetime(2017,04,01); 
G = [R RS C [S(t1); S(~t1)/2.0]];

G(G==0) = NaN;
N = sum(isfinite(G),2);
maxgas = medfiltSH(nanmax(G,[],2,'omitnan'),nDays);
medgas = medfiltSH(nanmedian(G,2),nDays);

figure('units','normalized','outerposition',[0 0 1/4 1]); 
ys = (2015:2021)';
lys = length(ys);
for i = 1:lys
    ax(i) = subplot(lys,1,i); 
    tI = t>=datetime(ys(i),01,01) & t< datetime(ys(i)+1,01,01); 
    ll = plot(ax(i),t(tI),medgas(tI),'linewidth',4); 
    ll.Color(4) = 0.5; zoom on; grid on;
    ax(i).YScale = 'log'; 
end
linkaxes(ax,'y');
title(ax(1),'"corrected" median gas: blue; "corrected" max gas: red');

figure(10); hold on; 
ys = (2015:2021)'; 
lys = length(ys); 
for i = 1:lys
    ax(i) = subplot(lys,1,i); hold on; tI = t>=datetime(ys(i),01,01) & t< datetime(ys(i)+1,01,01); 
    ll = plot(ax(i),t(tI),maxgas(tI),'linewidth',4); 
    ll.Color(4) = 0.5; zoom on; grid on; ax(i).YScale = 'log'; 
end
linkaxes(ax,'y');

figure(5);
axSmooth(2) = subplot(212);
plot(t,maxgas./medgas,'.'); 
zoom on; grid on; title('"corrected" Max / Median');

%%
close all; 
t1 = t <= datetime(2017,04,01); 
G = [R RS C [S(t1); S(~t1)/2.0]];
nDays = 30;
G(G==0) = NaN;
N = sum(isfinite(G),2);
maxgas = medfiltSH(nanmax(G,[],2,'omitnan'),nDays);
medgas = medfiltSH(nanmedian(G,2),nDays);

figure('units','normalized','outerposition',[0 0 1/4 1]);
hold on;
ys = (2015:2021)';
lys = length(ys);
for i = 1:lys
    subplot(211); hold on;
    tI = t>=datetime(ys(i),01,01) & t< datetime(ys(i)+1,01,01); 
    ll = plot(t(tI),medgas(tI),'linewidth',4); 
    ll.Color(4) = 0.5; zoom on; grid on;
end
ax = gca;
ax.YScale = 'log'; 

nDays = 7;
t1 = t <= datetime(2017,04,01); 
G = [R RS C [S(t1); S(~t1)/2.0]];
G(G==0) = NaN;
N = sum(isfinite(G),2);
maxgas = medfiltSH(nanmax(G,[],2,'omitnan'),nDays);
medgas = medfiltSH(nanmedian(G,2),nDays);

ys = (2015:2021)';
lys = length(ys);
for i = 1:lys
    subplot(212); hold on;
    tI = t>=datetime(ys(i),01,01) & t< datetime(ys(i)+1,01,01); 
    ll = plot(t(tI),medgas(tI),'linewidth',4); 
    ll.Color(4) = 0.5; zoom on; grid on;
end
ax = gca;
ax.YScale = 'log';
