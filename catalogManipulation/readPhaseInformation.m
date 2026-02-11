function [phaseStruct,n] = readPhaseInformation(event_file,GETP,originTime)
%
% [phaseStruct,n] = readPhaseInformation(event_file,phaseFlag,originTime)
%
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Saturday, Jul 27, 2019
%
% Rewrites on: 03FEB2026
%
%%
if nargin < 2
    GETP = true;
end

%% preallocate
fid = fopen(event_file);
phase_list = ["P";"Pn"];
if ~GETP
    phase_list = ["S";"Sn"];
end

%% get the job done
eofstat = feof(fid);
while ~eofstat
    buff = fgetl(fid);
    if buff == -1
        eofstat = true;
        continue;
    end
    tmp = sscanf(buff,'%d %s');
    charstr = char(tmp(2:end))';
    if ~strcmp(charstr,'Phase')
        if strcmp(charstr,'Station')
            %ive reached end of magnitude section
            break;
        end
        continue;
    end

    nphase = tmp(1);
    phaseStruct = populatePhaseStruct(nphase);
    t = NaT(nphase,1);
    sncls = repmat("",nphase,1);
    n = 0;
    phaseCount = 0;     % total number of all phases read
    buff = fgetl(fid);
    while phaseCount < nphase
        phaseCount = phaseCount + 1;
        C = textscan(buff,'%s %s %s %f %f %s %f %f %f %f %f %f %f %s %f %d %f %s');
        phase = string(C{6});
        lia = ismember(phase,phase_list);
        if ~lia
            buff = fgetl(fid);
            continue;
        end
        n = n + 1;
        stnm = C{1};
        stnm = string(stnm{1});
        ntwk = C{2};
        ntwk = string(ntwk{1});
        chan = C{3};
        chan = chan{1};
        chan = string(chan(1:end));

        t_ = datetime(C{7},C{8},C{9},C{10},C{11},C{12});
        tt = milliseconds(t_ - originTime);
        residual = C{13};
        weight = C{15};
        polarity = C{16};
        takeOffAngle = C{17};
        khole = C{18};
        if ~isempty(khole)
            khole = khole{1};
        else
            khole = "";
        end
        phaseStruct(n).stnm = stnm;
        phaseStruct(n).ntwk = ntwk;
        phaseStruct(n).chan = chan;
        phaseStruct(n).locid = khole;
        phaseStruct(n).dist = C{4};
        phaseStruct(n).azimuth = C{5};
        phaseStruct(n).t = t_;
        phaseStruct(n).res = residual;
        phaseStruct(n).wt = weight;
        phaseStruct(n).polarity = polarity;
        phaseStruct(n).tt = tt;
        phaseStruct(n).toa = takeOffAngle;
        t(n) = t_;
        sncls(n) = sprintf("%s.%s.%s.%s",stnm,ntwk,chan,khole);
        buff = fgetl(fid);
    end
    break;
end
fclose(fid);

%%
phaseStruct = phaseStruct(1:n);
if ~n
    return;
end
sncls = sncls(1:n);
t = t(1:n);
[~,uI] = unique(sncls);         % remove dups
phaseStruct = phaseStruct(uI);  % remove dups
t = t(uI);                      % remove dups
[~,sI] = sort(t);               
phaseStruct = phaseStruct(sI);
n = length(phaseStruct);