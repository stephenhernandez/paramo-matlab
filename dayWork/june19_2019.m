% %%june19_2019
clear; close all; clc;
% lbc = LibComCat();
% T = (1800:2020)';
% lw = length(T);
% t = [];
% eqlat = [];
% eqlon = [];
% eqdepth = [];
% eqmag = [];
% magType = [];
% status = [];
% code = [];
% 
% %%
% n = 1;
% while n <= lw-1
%     disp(T(n));
%     events = lbc.getEventData('starttime',datenum(T(n),01,01),'endtime',datenum(T(n+1),01,01) ...);
%         ,'xmin',-180,'xmax',180,'ymin',-90,'ymax',90,'minmag',5,'maxmag',10);
%     U = libcomcat2struct(events);
%     t_ = pull(U,'t');                        t = [t; t_];
%     eqlat_ = pull(U,'lat');                  eqlat = [eqlat; eqlat_];
%     eqlon_ = pull(U,'lon');                  eqlon = [eqlon; eqlon_];
%     eqmag_ = pull(U,'mag');                  eqmag = [eqmag; eqmag_];
%     eqdepth_ = pull(U,'depth');              eqdepth = [eqdepth; eqdepth_];
%     magType_ = string(pull(U,'magType'));    magType = [magType; magType_];
%     status_ = string(pull(U,'status'));      status = [status; status_];
%     code_ = string(pull(U,'code'));          code = [code; code_];
%     n = n+1;
% end
% 
% %%
% [t,ia] = unique(t);
% eqlat = eqlat(ia);
% eqlon = eqlon(ia);
% eqdepth = eqdepth(ia);
% eqmag = eqmag(ia);
% magType = magType(ia);
% status = status(ia);
% code = code(ia);

%%
cd ~/igdata/
updateGlobalCatalog();
load globalCatalog.mat
ia = t >= datetime(1990,01,01) & eqmag >= 5;
t = t(ia);
eqlat = eqlat(ia);
eqlon = eqlon(ia);
eqdepth = eqdepth(ia);
eqmag = eqmag(ia);
magType = magType(ia);
status = status(ia);
code = code(ia);

%%
figure();
S = scatter(eqlon,eqlat,0.25*exp(eqmag),log10(eqdepth),'filled');
colorbar;
axis equal;
zoom on;
S.MarkerFaceAlpha = 0.5;
caxis([1 3]);

refEllipse = referenceEllipsoid('wgs84');
[stla,stlo,stel] = metaDataFromStationList("VCH1");
d = distance(eqlat,eqlon,stla,stlo,refEllipse)*1e-3;
degs = km2deg(d);
Dmin = 200;
dI = find(true(size(d)));
d2 = degs(dI);
log10Amp = eqmag(dI)-(0.5*log10(sind(degs(dI))))-...
    (0.0031*degs(dI))+log10(0.6./(20*sqrt(degs(dI))))+0.43 -0; %-0 for nm, -3 for microns, -6 to mm, -9 for meters

phasevelocity=3.5e3*1e9; % c in nm
predAmpRussell = 10.^log10Amp;
russellStrain = (2*pi*predAmpRussell) / (20*phasevelocity); %amplitude is now strain

largestN = 400;
magFact = 0.25;
alphaValue = 0.5;
[tmpampSort,aI] = sort(russellStrain,'descend');

%%
figure();
S = scatter(eqlon(dI(aI(1:largestN))),eqlat(dI(aI(1:largestN))),magFact*exp(eqmag(dI(aI(1:largestN)))),datenum(t(dI(aI(1:largestN)))),'filled');
colorbar;
S.MarkerFaceAlpha = alphaValue;
S.MarkerEdgeColor = 'k';
hold on;
plot(stlo+0,stla,'^','markersize',20);
zoom on;

%%
figure();
plot(t(dI(aI(1:largestN))),d(dI(aI(1:largestN))),'o');
zoom on;

%%
figure();
semilogy(t,russellStrain,'.');
hold on;
semilogy(t(dI(aI(1:largestN))),russellStrain(dI(aI(1:largestN))),'o'); zoom on;

%%
figure();
plot(t,1:length(t),'o');
zoom on;

%%
figure();
plot(t,eqmag,'.');
hold on;
plot(t(dI(aI(1:largestN))),eqmag(dI(aI(1:largestN))),'o');
zoom on;

%%
figure();
plot(d,russellStrain,'k.');
hold on;
S = scatter(d(dI(aI(1:largestN))),russellStrain(dI(aI(1:largestN))),magFact*exp(eqmag(dI(aI(1:largestN)))),datenum(t(dI(aI(1:largestN)))),'filled');
colorbar;
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
zoom on;
S.MarkerFaceAlpha = alphaValue;
