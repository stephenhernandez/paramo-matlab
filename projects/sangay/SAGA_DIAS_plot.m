clear; close all;
newFs = 1;
npoles = 1;
lfc = 1/320;
hfc = 1/5;

S = loadWaveforms(datetime(2021,04,29),1,"SAGA","HHN");
Scut = cutWaveforms(S,dateshift(S(1).ref,'start','day')+hours(12)+minutes(00),0,hours(2));
Scut = detrendWaveforms(Scut);
Scut(2) = filterWaveforms(Scut(1),-inf,hfc,1);
Scut(3) = scaleWaveforms(transferWaveforms(resampleWaveforms(Scut(1),newFs),lfc,hfc,npoles,newFs,'disp'),1e6);
Scut(4) = scaleWaveforms(transferWaveforms(resampleWaveforms(Scut(1),1),lfc,hfc,npoles,newFs,'acc'),1e9/9.81);

%%
close all;
[~,ax] = plotWaveforms(Scut);
titles = ["raw";"raw low passed at 0.1 hz. (pass = 1, poles = 1)";...
    "displacement (microns); Filter: 20 - 320 sec., 1 pass, 2 poles";...
    "tilt (nanradians); Filter: 20 - 320 sec., 1 pass, 2 poles"];

for i = 1:length(ax)
    disp(i);
    ax(i).Title.String = titles(i);
end