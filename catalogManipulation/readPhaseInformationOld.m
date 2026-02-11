function [origmag,nPphases,Pphases] = readPhaseInformationOld(event_file,origyyyy,origmm,origdd,orighh,origmmm,origsec)
fid = fopen(event_file);
eofstat = feof(fid);

%preallocate
origmag = NaN;
magCount = 0;
phaseCount = 0;
nFlag = false;
pFlag = false;
origin_time = datetime(origyyyy,origmm,origdd,orighh,origmmm,origsec);
n = 0;

%get the job done
while ~eofstat
    buff = fgetl(fid);
    if buff ~= -1
        tmp = sscanf(buff,'%d %s');
        charstr = char(tmp(2:end))';
        if strcmp(charstr,'Network')
            nmags = tmp(1);
            nFlag = true;
        elseif strcmp(charstr,'Phase')
            nphase = tmp(1);
            Pphases = populatePhaseStruct(nphase);
            pFlag = true;
            nFlag = false;
        elseif strcmp(charstr,'Station')
            pFlag = false;
        elseif nFlag
            if magCount <= nmags
                C = textscan(buff,'%s %f %d %s');
                if strcmp(C{4},'preferred')
                    origmag = C{2};
                end
                magCount = magCount + 1;
            end
        elseif pFlag
            phaseCount = phaseCount + 1;
            if phaseCount > 1 && phaseCount < nphase + 2 && nphase > 0
                C = textscan(buff,'%s %s %f %d %s %f:%f:%f %f %s %f %s');
                phase = C{5};
                if strcmp(phase,'P')
                    n = n + 1;
                    stnm = C{1}; Pphases(n).stnm = string(stnm{1});
                    ntwk = C{2}; Pphases(n).ntwk = string(ntwk{1});
                    Pphases(n).dist = C{3};
                    Pphases(n).azimuth = C{4};
                    hh = C{6};
                    mmm = C{7};
                    sec = C{8};
                    residual = C{9};
                    weight = C{11};
                    if isempty(residual)
                        residual = NaN;
                        weight = NaN;
                    end
                    %Pphases(n).hh = hh;
                    %Pphases(n).mmm = mmm;
                    %Pphases(n).sec = sec;
                    Pphases(n).res = residual;
                    Pphases(n).wt = weight;
                    first_try = datetime(origyyyy,origmm,origdd,hh,mmm,sec);
                    theoretical_travel_time = seconds(first_try - origin_time);
                    if theoretical_travel_time < 120
                        Pphases(n).t = first_try;
                    else
                        doy = ymd2doy(origyyyy,origmm,origdd);
                        doy = doy+1;
                        [newYear,newMonth,newDay] = ydoy2ymd(origyyyy,doy);
                        Pphases(n).t = datetime(newYear,newMonth,newDay,hh,mmm,sec);
                    end
                end
            end
        end
    end
    eofstat = feof(fid);
end
fclose(fid);

%% now sort
Pphases = Pphases(1:n);
t = pull(Pphases,'t');
[~,tI] = sort(t);
Pphases = Pphases(tI);
nPphases = n;

%% these bits are weird, possibly unnecessary s.h., sat. 03 aug., 2019
% if nphase
%     [status,Pphases] = removePphaseDuplicates(Pphases)
%     while status
%         [status,Pphases] = removePphaseDuplicates(Pphases)
%     end
% end
% nPphases = min([nphase length(Pphases)]);