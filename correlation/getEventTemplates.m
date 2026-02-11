function status = getEventTemplates(P,masterStructure,lead,lag)

if nargin < 4
    lag = 2.5;
end
if nargin < 3
    lead = -0.5;
end
status = 0;
nEvents = length(P);

% Load master structure
existFlag = exist([masterStructure,'.mat'],'file');
if existFlag
    disp('Loading master structure');
    load(masterStructure);
else
    disp('Master data structure does not exist: try again');
    status = 1;
    return
end

% Proceed to making templates per event
Fs = 1/S(1).delta; %Since S is synched, Fs is the same for all streams...
winLen = Fs*(lag-lead); %Total windlow length
for i = 1:nEvents
    templateStructure = P(i).id;
    origin = P(i).OabsTime;
    magnitude = P(i).origmag;
    disp(['Working on: ',templateStructure])
    Pphase = P(i).Pphases;
    [Pphase,nPhases] = reducePhases(Pphase);
    
    if nPhases > 0
        % Pre-Allocate
        T = populateTemplateStructure;
        T(1).mintt = Pphase(1).PabsTime + 1;
        tmpShifts = findShifts(Pphase,Fs);
        currentTimeVector = S(1).ref + (0:S(1).npts-1)*S(1).delta/86400; % Since S is synched, we only have to do this once
        startTime = currentTimeVector(1);
        endTime = currentTimeVector(end);
        endTime = endTime - (lag/86400);
        
        % Get templates into matrix d
        n = 0;
        for j =  1:length(S)
            currentSname = S(j).kstnm;
            currentScmpnm = S(j).kcmpnm;
            for k =1:nPhases
                Pphasename = Pphase(k).stnm;
                Pphasetime = Pphase(k).PabsTime + (lead/86400);
                %check if any amount of data from this station exist
                if strcmp(currentSname,Pphasename)
                    %disp(['Large stream begins at: ',datestr(startTime)])
                    %disp(['Large stream ends at: ',datestr(endTime)])
                    %disp(['P time is: ', datestr(Pphasetime)])
                    % Check that data can actually be recovered
                    if Pphasetime >= startTime && Pphasetime <= endTime
                        n = n + 1;
                        dtmp = S(j).d;
                        startIndex = find(Pphasetime-currentTimeVector <= 0,1,'first');
                        %disp(currentSname)
                        %disp(['Start Index: ',num2str(startIndex)])
                        T(n).d = dtmp(startIndex:startIndex+winLen-1);
                        T(n).shift = tmpShifts(k);
                        T(n).chan = currentScmpnm;
                        T(n).stnm = currentSname;
                        if tmpShifts(k) == 0
                            mintt = Pphasetime;
                        end
                        T(n).origin = origin;
                        T(n).id = templateStructure;
                        T(n).mag = magnitude;
                    else
                        disp(strcat('Data do not exist for: ',currentSname,currentScmpnm));
                        continue;
                    end
                end
            end
        end
        
        if n > 0
            disp(['Saving Event ID: ',templateStructure]);
            for ii = 1:n
                T(ii).mintt = mintt;
            end
            save(templateStructure,'T','-v7.3');
            clear T d shifts
        else
            disp(['Event ID ',templateStructure,' doesnt have data. Skipping...'])
        end
    else
        disp(['Event ID ',templateStructure,' doesnt have phases. Skipping...']);
        continue;
    end
end