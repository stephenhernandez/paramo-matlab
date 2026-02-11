function [origmag,nSphases,Sphases] = readSphaseInformationOld(event_file,origyyyy,origmm,origdd,orighh,origmmm,origsec)
fid = fopen(event_file);
eofstat = feof(fid);

%preallocate
origmag = NaN;
magCount = 0;
phaseCount = 0;
nFlag = false;
sFlag = false;
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
            nphase = tmp(1)+1;
            Sphases = populatePhaseStruct(nphase);
            sFlag = true;
            nFlag = false;
        elseif strcmp(charstr,'Station')
            sFlag = false;
        elseif nFlag
            if magCount <= nmags
                C = textscan(buff,'%s %f %d %s');
                if strcmp(C{4},'preferred')
                    origmag = C{2};
                end
                magCount = magCount + 1;
            end
        elseif sFlag
            phaseCount = phaseCount + 1;
            if phaseCount > 1 && phaseCount < nphase + 1
                C = textscan(buff,'%s %s %f %d %s %f:%f:%f %f %s %f %s');
                phase = C{5};
                if strcmp(phase,'S')
                    n = n + 1;
                    stnm = C{1}; Sphases(n).stnm = stnm{1};
                    ntwk = C{2}; Sphases(n).ntwk = ntwk{1};
                    Sphases(n).dist = C{3};
                    Sphases(n).azimuth = C{4};
                    hh = C{6};
                    mmm = C{7};
                    sec = C{8};
                    Sphases(n).hh = hh;
                    Sphases(n).mmm = mmm;
                    Sphases(n).sec = sec;
                    Sphases(n).res = C{9};
                    Sphases(n).wt = C{11};
                    first_try = datetime(origyyyy,origmm,origdd,hh,mmm,sec);
                    theoretical_travel_time = seconds(first_try - origin_time);
                    if theoretical_travel_time < 120
                        Sphases(n).t = first_try;
                    else
                        doy = ymd2doy(origyyyy,origmm,origdd);
                        doy = doy+1;
                        [newYear,newMonth,newDay] = ydoy2ymd(origyyyy,doy);
                        Sphases(n).t = datetime(newYear,newMonth,newDay,hh,mmm,sec);
                    end
                else
                    continue
                end
            end
        end
    end
    eofstat = feof(fid);
end
fclose(fid);

%%
nSphases = n;
Sphases = Sphases(1:n);

%%
% [status,Sphases,nSphases] = removeSphaseDuplicates(Sphases);
% while status
%     [status,Sphases,nSphases] = removeSphaseDuplicates(Sphases);
% end