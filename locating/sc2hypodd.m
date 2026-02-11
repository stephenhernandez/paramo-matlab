function sc2hypodd(E,outfile,depCorr,minUsuablePhases,minSPhases)
%
% sc2hypodd(E,outfile,depCorr,minUsuablePhases,minSPhases)
%

%
% Rewrites on 02FEB2026
%

%%
if nargin < 3
    depCorr = 0;
end

if nargin < 4
    minUsuablePhases = 3;
end

if nargin < 5
    minSPhases = 0;
end
%%
nS = pull(E,"nSphases");
nP = pull(E,"nPphases");
phasesTot = nS + nP;
npI = phasesTot > minUsuablePhases & nS >= minSPhases;
E = E(npI);

%%
t = pull(E,"t");
lat = pull(E,"lat");
lon = pull(E,"lon");
depth = pull(E,"depth");
laterr = pull(E,"laterr");
lonerr = pull(E,"lonerr");
mag = pull(E,"mag");
rms = pull(E,"rms");
idOrig = pull(E,'id',"");

%%
deptherr = (laterr + lonerr)/2;
horizErr = 1*sqrt(laterr.^2 + lonerr.^2);

%%
formatSpec1 = '# %d %02d %02d %02d %02d %5.2f %8.4f %8.4f %8.2f %5.2f %5.2f %6.2f %5.2f %s %s';
formatSpec2 = '%s %f %f %1s';
N = length(t);
strMain = [];
pStnms = [];
phase_arrival_time = [];
phase_residual = [];
phaseStr = [];

for i = 1:N
    thisID = idOrig(i);
    if horizErr(i) > 1000
        fprintf("cant process %s, continuing\n",thisID);
        continue;
    end
    origtime = t(i);
    [yyyy_,month_,day_,hour_,minute_,seconds_] = datevec(t(i));

    idStr = sprintf("%d%05d",yyyy_,i);
    str1 = compose(formatSpec1,yyyy_,month_,day_,hour_,minute_,seconds_,...
        lat(i),lon(i),depth(i)+depCorr,mag(i),horizErr(i),deptherr(i),...
        rms(i),idStr,thisID);

    E_ = E(i);
    nP_ = E_.nPphases;
    if ~nP_
        fprintf("cant process %s, not enough P phases\n",thisID);
        continue;
    end

    Pphases = E_.Pphases;
    phaseStr = repmat('P',nP_,1);
    pStnms = string(pull(Pphases,"stnm"));
    phase_arrival_time = pull(Pphases,"t");
    phase_residual = pull(Pphases,"res");

    nS_ = E_.nSphases;
    if nS_
        Sphases = E_.Sphases;
        phaseStr = [phaseStr; repmat("S",nS_,1)];
        pStnms = [pStnms; string(pull(Sphases,"stnm"))];
        phase_arrival_time = [phase_arrival_time; pull(Sphases,"t")];
        swts = pull(Sphases,"res");
        phase_residual = [phase_residual;swts];
    end

    phase_arrival_time = seconds(phase_arrival_time - origtime);
    [phase_arrival_time,sI] = sort(phase_arrival_time);
    phaseStr = phaseStr(sI);
    pStnms = pStnms(sI);
    phase_residual = phase_residual(sI);
    %1-(abs(zscore(phase_residual))/4)
    new_weight1 = floor(-phase_residual + abs(phase_residual)/mad(phase_residual,0)); %use 0 method for MAD... less strict and ultimately better
    new_weight = 1-(new_weight1/4);
    %table(pStnms,phase_arrival_time,new_weight1,new_weight,phaseStr)
    %new_weight_alternative = 1-(abs(zscore(phase_residual))/4); %alternative definition

    %
    %lia = ismember(pStnms,["T01";"T02";"T03";"T04";"T05";"T06";"T07";"T08";"T09";"T10";"T11";"T12"]);
    sI = new_weight >= 0; % & ~lia;
    if sum(sI) < minUsuablePhases
        continue;
    end
    phase_arrival_time = phase_arrival_time(sI);
    new_weight = new_weight(sI);
    phaseStr = phaseStr(sI);
    pStnms = pStnms(sI);

    str2 = compose(formatSpec2,pStnms,phase_arrival_time,new_weight,phaseStr);
    strMain = [strMain; string(str1); string(str2)];
    % test_cmd = sprintf("\\cp ~/phaseInformationSC5/*/%s* ~/masa/relocation/for_abigail/\n",thisID);
    % unix(test_cmd)
end

%%
disp(outfile)
fileID = fopen(outfile,'w');
fprintf(fileID,'%s\n',strMain);
fclose(fileID);