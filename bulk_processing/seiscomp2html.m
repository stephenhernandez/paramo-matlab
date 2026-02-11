% function S = seiscomp2html(fname)
clear; close all; clc;
%fname = '~/phaseInformationSC3/2019/igepn2019clmg.txt';
fname = 'igepn2013upau'; %igepn2019ixzp.txt';

%%
PWD = pwd;
E = readSCBulletin(fname); %'~/phaseInformationSC3/2016/igepn2016hnmu.txt');
nPphases = E.nPphases;
t = E.t;
eqlat = E.lat;
eqlon = E.lon;
eqmag = E.mag;
if eqmag == -999
    eqmag = 7.8;
end
startTime = -5;
refEllipse = referenceEllipsoid('wgs84');
beginTime = t + seconds(startTime);

%%
n = 0;
S = populateWaveforms(nPphases);
endTime = NaT(nPphases,1);
d = NaN(nPphases);

if nPphases
    pphases = E.Pphases;
    for i = 1:nPphases
        pt = pphases(i).t;
        stnm = pphases(i).stnm;
        charchan = char(pphases(i).chan);
        netwk = pphases(i).ntwk;
        [stla,stlo] = metaDataFromStationList(string(stnm));
        d_ = distance(eqlat,eqlon,stla,stlo,refEllipse)*1e-3;
        logt = -1.04 + 0.44*eqmag + 0.19*log(d_); %these values are from Hernandez and Cotton (2000)
        dur_ = 12*exp(logt); %105;
        
        S_ = extractWaveforms(beginTime,pt+seconds(dur_),stnm,[charchan(1:2),'Z'],netwk,"",false,0);
        if ~isnat(S_.ref)
            n = n+1;
            S(n) = S_;
            endTime(n) = pt+seconds(dur_);
            d(n) = d_;
        end
    end
end

%%
cd ~/events/html/
id = E.id;

if ~exist(id,'dir')
    mkdir(id)
end
cd(id);

S = S(1:n);
S = intWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),0.5,8),0.01));
fprintf('\n\n%s\n\n',id);

close all;
ax = plotWaveforms(S(1));
h = gcf;
h.Visible = 'off';
set(h,'defaultLegendAutoUpdate','off');
for i = 1:n
    S_ = S(i);
    stnm_ = char(S_.kstnm);
    net_ = char(S_.knetwk);
    chan_ = char(S_.kcmpnm);
    locID_ = char(S_.khole);
    cla(ax);
    [ax,hp] = plotWaveforms(ax,S_);
    hp.LineWidth = 1;
    hold(ax,'on');
    plot([t t],[S_.depmin S_.depmax],'k-','linewidth',1);
    printName = [id,'_',stnm_,'_',net_,'_',locID_,'_',chan_];
    fprintf('%s.%s.%s.%s: %s - %s, %5.2f km\n',stnm_,net_,locID_,chan_,datestr(beginTime),datestr(endTime(i)),d(i));
    print(printName,'-dpng','-r100');
end
cd(PWD);
