clear; close all;
cd ~/data/payg/;

dayStart = datetime(2023,01,350);
dayEnd = datetime(2024,01,040); %datetime('now');
dayInc = 1;
%dayVec = (dayStart:dayInc:dayEnd)';
dayVec = (dayEnd:-dayInc:dayStart)';
ldays = length(dayVec);

%%
cmpnmList = ["BHZ";"BH1";"BH2"];
lCmpnms = length(cmpnmList);
for i = 1:ldays
    day_ = dayVec(i);
    disp(day_);
    tstart = datestr(day_,0);
    tend = datestr(day_+1,0);
    for j = 1:lCmpnms
        cmpnm_ = char(cmpnmList(j));
        S = irisFetch.Traces('IU','PAYG','00',cmpnm_,tstart,tend);
        if isempty(S)
            continue;
        end

        S = iris2sh(S);
        S = S(:);
        try
            S = mergeWaveforms(S);
        catch
            fprintf('some issue with merging, moving on...\n');
            continue;
        end

        if ~isscalar(S)
            continue;
        end

        S.d = single(S.d);

        ref = S.ref;
        yyyy = year(ref);
        doy = day(ref,'doy');

        fname = sprintf("%s.%s.%s.%s.D.%d.%03d",S.knetwk,S.kstnm,S.khole,S.kcmpnm,yyyy,doy);
        mkmseed(char(fname),S.d,datenum(ref),1./S.delta,11,4096*2);
        fprintf('\n');
        %S = struct2table(S,'AsArray',true);
        %save(fname,'S','-v7.3');
    end
end
