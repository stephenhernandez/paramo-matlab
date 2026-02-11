function scthree2hypodd(E,outfile,depCorr,minUsuablePhases)
if nargin < 3; depCorr = 0; end
if nargin < 4; minUsuablePhases = 3; end
NP = pull(E,'nSphases') + pull(E,'nPphases');
npI = NP > minUsuablePhases;
E = E(npI);

%%
t = pull(E,'t');
lat = pull(E,'lat');
lon = pull(E,'lon');
depth = pull(E,'depth');
laterr = pull(E,'laterr');
lonerr = pull(E,'lonerr');
mag = pull(E,'mag');
rms = pull(E,'rms');

%%
deptherr = (laterr + lonerr)/2;
horizErr = 1*sqrt(laterr.^2 + lonerr.^2);

%%
formatSpec1 = '%s %d %02d %02d %02d %02d %5.2f %8.4f %8.4f %8.2f %5.2f %5.2f %6.2f %5.2f %s';
formatSpec2 = '%s %11.3f %3d %7s';
N = length(t);
strMaster = [];
pStnms = [];
ptt = [];
wt = [];
phaseStr = [];

for i = 1:N
    if horizErr(i) > 1000
        disp(['cant process',num2str(i)]);
        continue;
    end
    origtime = t(i);
    [yyyy_,month_,day_,hour_,minute_,seconds_] = datevec(t(i));

    idStr = num2str(yyyy_);
    if i < 10
        idStr = [idStr,'0000',num2str(i)];
    elseif i < 100
        idStr = [idStr,'000',num2str(i)];
    elseif i < 1000
        idStr = [idStr,'00',num2str(i)];
    elseif i < 10000
        idStr = [idStr,'0',num2str(i)];
    else
        idStr = [idStr,num2str(i)];
    end

    str1 = compose(formatSpec1,'# ',...
        yyyy_,month_,day_,hour_,minute_,seconds_,lat(i),lon(i),depth(i)+depCorr,mag(i),horizErr(i),deptherr(i),rms(i),idStr);
    if E(i).nPphases
        Pphases = E(i).Pphases;
        phaseStr = repmat('P',E(i).nPphases,1);
        pStnms = string(pull(Pphases,'stnm'));
        ptt = pull(Pphases,'t');
        wt = pull(Pphases,'wt');
    end
    if E(i).nPphases && E(i).nSphases
        Sphases = E(i).Sphases;
        phaseStr = [phaseStr;repmat('S',E(i).nSphases,1)];
        pStnms = [pStnms;string(pull(Sphases,'stnm'))];
        ptt = [ptt;pull(Sphases,'t')];
        wt = [wt;pull(Sphases,'wt')];
    end
    ptt = seconds(ptt - origtime);
    [ptt,sI] = sort(ptt);
    phaseStr = phaseStr(sI);
    pStnms = pStnms(sI);
    wt = wt(sI);
    str2 = compose(formatSpec2,pStnms,ptt,wt,phaseStr);
    strMaster = [strMaster;string(str1);string(str2)];
end

%%
fileID = fopen(outfile,'w');
fprintf(fileID,'%s\n',strMaster);
fclose(fileID);
