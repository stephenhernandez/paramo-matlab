clear; close all; clc; 
[eqType,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,nPhases,nMLv,...
timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
locMethod,earthModel,creationTime] = readCat1();
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(eqlat,eqlon,-0.53,-78.49,refEllipse)*1e-3;
dI = d_ <= 5;

figure(); plot(t(dI),1:sum(dI),'.'); zoom on;
figure(); ss = scatter(eqlon(dI),eqlat(dI),3*exp(eqmag(dI)),datenum(t(dI)),'filled'); zoom on; grid on; c = colorbar;
axis equal;
ss.MarkerFaceAlpha = 0.5;
ss.MarkerEdgeColor = 'k';
ss.MarkerEdgeAlpha = 0.25;
c.TickLabels = datestr(c.Ticks);


SZ = extractWaveforms(t(dI)+seconds(10),seconds(30),"BTER","HHZ","EC","",true,true,1);
SN = extractWaveforms(t(dI)+seconds(10),seconds(30),"BTER","HHN","EC","",true,true,1);
SE = extractWaveforms(t(dI)+seconds(10),seconds(30),"BTER","HHE","EC","",true,true,1);


% SZ = extractWaveforms(t(dI),seconds(20),"TOMA","BHZ","EC","",true,true,1);
% SN = extractWaveforms(t(dI),seconds(20),"TOMA","BHN","EC","",true,true,1);
% SE = extractWaveforms(t(dI),seconds(20),"TOMA","BHE","EC","",true,true,1);

%%
close all; 
rIZ = isnat(pull(SZ,'ref')) | ~isfinite(pull(SZ,'depmen'));
rIN = isnat(pull(SN,'ref')) | ~isfinite(pull(SN,'depmen'));
rIE = isnat(pull(SE,'ref')) | ~isfinite(pull(SE,'depmen'));

rI = rIZ | rIN | rIE;
SZ(rI) = [];
SN(rI) = [];
SE(rI) = [];

lfc = 4;
hfc = 8;
tw = 0.04;
dZ = double(pull(intWaveforms(filterWaveforms(taperWaveforms(differentiateWaveforms(SZ),tw),lfc,hfc))));
dN = double(pull(intWaveforms(filterWaveforms(taperWaveforms(differentiateWaveforms(SN),tw),lfc,hfc))));
dE = double(pull(intWaveforms(filterWaveforms(taperWaveforms(differentiateWaveforms(SE),tw),lfc,hfc))));

d = normalizeWaveforms([normalizeWaveforms(dZ); normalizeWaveforms(dN); normalizeWaveforms(dE)]);
[maxccp,plags,maxccn,nlags] = doccFreqCircShift(d(:,peak2rms(d)' >= 4),true);

d2 = apply_vdcc(d); 
tmp_stack = plot_family(d2,(1:size(d2,2))',5,100);

[newEvents,newFamilies] = pruneAndMergeEvents(d,maxccp,plags,[0.8 0.7 0.7],'weighted');
d2 = newEvents;
tmp_stack = plot_family(d2,(1:size(d2,2))',5,100);

%newEvents = apply_vdcc(newEvents);
d2 = newEvents; 
tmp_stack = plot_family(d2,(1:size(d2,2))',5,100);


% figure(); imagesc(squareform(maxccp)); colorbar; zoom on;
% 
% d2 = newEvents; 
% tmp_stack = plot_family(d2,(1:size(d2,2))',5,100);
% 
[U,SVs,V] = svd(newEvents,0);
d2 = U; tmp_stack = plot_family(d2,(1:size(d2,2))',5,100);
figure(); semilogy(diag(SVs),'.'); zoom on; grid on;

%%
U_ = U;
clearvars -except U_;
load('~/research/now/machachi/bter_svd_basis_functions');
clear U
U(:,1,:) = U_(1:3000,:); % U_(2001:4000,:) U_(4001:end,:)];
U(:,2,:) = U_(3001:6000,:); % U_(2001:4000,:) U_(4001:end,:)];
U(:,3,:) = U_(6001:end,:); % U_(2001:4000,:) U_(4001:end,:)];
clear U_;
save('~/research/now/machachi/bter_svd_basis_functions');

