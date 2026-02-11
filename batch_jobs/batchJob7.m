clear; close all;
cd ~/research/now/cotopaxi/
load BREF_templates_3.mat
T2 = table2struct(T2);
T = reshape(T2,sizeT);
clear T2

lfc = 3/8;
hfc = 3;
newFs = 50;

dayStart = datetime(2022,11,24);
dayEnd = datetime(2022,11,24);

dayVec = (dayEnd:-1:dayStart)';
lDays = length(dayVec);

for i = 1:lDays
    day_ = dayVec(i);
    C = loadWaveforms(day_,1,["BREF";"BTAM"],["BHZ";"BHN";"BHE"],"EC");
    if isnat(C(1).ref)
        continue;
    end

    C = cutWaveforms(C,dateshift(C(1).ref,'start','day')+hours(1)+minutes(20)+seconds(0),0,minutes(10));

    Cf = syncWaveforms(filterWaveforms(detrendWaveforms(C),lfc,hfc));
    Cf = resampleWaveforms(Cf,newFs);
    Cf = nanGapWaveforms(Cf,0);

    [ccMain,tMain,ampMain,dMag,evidMain,magMain,templateNumber,madMain] = ...
        templateSearch(Cf,T(:,753),false,true,true);
end

% %C = loadWaveforms(datetime(2021,11,27),2,["BREF";"BTAM"],["BHZ";"BHN";"BHE"],"EC");
% close all; [indiv_events,tabs,NCC,z2p,Neff,p2rms,kurt,ccnorm,t] = ...
% subspaceDetector(C,0.2,basisFunctionName,15,...
% 23,1e4,20,false,true,true,true,false);
% %