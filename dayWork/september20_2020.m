clear; close all; clc;

%%
S = loadWaveforms(datetime(2020,06,08),2,"PUYO","HHZ");
S(2) = loadWaveforms(datetime(2020,09,20),1,"PUYO","HHZ");

%%
Sf = detrendWaveforms(filterWaveforms(S,0.6,1.2));
Senv1 = envelopeWaveforms(Sf);

%%
plotWaveforms(Sf,[],[],[],[],true);
[~,ax] = plotWaveforms(Senv1,[],[],[],[],true);
ax(1).YScale = 'log'; ax(2).YScale = 'log';

%%
Sf(1).d = (Sf(1).d).^2;
Sf(2).d = (Sf(2).d).^2;

%%
Senv2 = Senv1;
Senv2(1).d = (Senv2(1).d).^2;
Senv2(2).d = (Senv2(2).d).^2;


%%
[~,ax] = plotWaveforms(Senv2,[],[],[],[],true);
ax(1).YScale = 'log'; ax(2).YScale = 'log';

%%
Senv3 = Senv2;
Senv3(1).d = zpkFilter(Senv3(1).d,-inf,1/100,1,1,1);
Senv3(2).d = zpkFilter(Senv3(2).d,-inf,1/100,1,1,1);

plotWaveforms(Senv3,[],[],[],[],true);

energy09jun = cutWaveforms(Senv3(1),dateshift(Senv3(1).ref,'start','day')+hours(24)+minutes(00)+seconds(00),0,hours(5),true);
energy20sep = cutWaveforms(Senv3(2),dateshift(Senv3(2).ref,'start','day')+hours(07)+minutes(30)+seconds(00),0,hours(5),true);


%%
[~,ax] = plotWaveforms([energy09jun;energy20sep],[],[],[],[],true); ax(1).YScale = 'log'; ax(2).YScale = 'log';
disp(sum(energy20sep.d)./sum(energy09jun.d))

%%
%ax(1).YScale = 'log'; %ax(2).YScale = 'log';
%[~,ax] = plotWaveforms(Senv3(2));
%ax(1).YScale = 'log'; %ax(2).YScale = 'log';