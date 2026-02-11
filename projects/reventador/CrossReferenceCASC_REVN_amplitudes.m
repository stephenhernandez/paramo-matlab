clear; close all;

cd ~/research/now/reventador/
load LizGoodData_03DEC2021.mat

figure(); semilogy(tGood,aGood,'.'); zoom on;
tI = tGood >= datetime(2018,01,01) & tGood < datetime(2019,01,01);
tGood = tGood(tI);
aGood = aGood(tI);

%%
lfc = 1/4;
hfc = 2;

revn_amps = NaN(size(tGood));
revs_amps = revn_amps;
casc_amps = revn_amps;

revn_t = NaT(size(tGood));
revs_t = revn_t;
casc_t = revn_t;

dayStart = datetime(2018,01,01);
dayEnd = datetime(2018,12,31);
days = (dayStart:dayEnd)';
lDays = length(days);

%%
si = 1;
for i = 1:lDays-1
    day_ = days(i);
    disp(day_);
    tI = tGood >= day_ & tGood < day_+1;
    sumti = sum(tI);
    if ~sumti
        fprintf('no events on this day\n');
        continue;
    end

    %
    ei = si+sumti-1;
    SCASC = loadWaveforms(day_,1,"CASC","HHZ","EC");
    if isnat(SCASC.ref)
        continue;
    end

    SREV = loadWaveforms(day_,1,["REVN";"REVS"],"HHZ","EC");
    if isnat(SREV(1).ref)
        continue;
    end

    tToday = tGood(tI);

    Sf = transferWaveforms(detrendWaveforms(differentiateWaveforms(SCASC)),lfc,hfc,4,25,'disp',true,false);
    Scut = cutWaveforms(Sf,tToday,0,90);

    if isnat(Scut(1).ref)
        disp('something went wrong');
        continue;
    end

    %
    try
        ampsTmp = 0.5*peak2peak(double(pull(Scut)))';
        refsTmp = pull(Scut,'ref');
        lcut = length(ampsTmp);

        casc_amps(si:si+lcut-1) = ampsTmp;
        casc_t(si:si+lcut-1) = refsTmp;
        whos
    catch
        fprintf('couldnt strip casc of its amplitudes, continuing\n');
        continue;
    end

    %
    lRev = length(SREV);
    Sf = transferWaveforms(detrendWaveforms(differentiateWaveforms(SREV)),0.75,3,4,25,'disp',true,false);
    for j = 1:lRev
        Scut = cutWaveforms(Sf(j),tToday,0,90);
        try
            ampsTmp = 0.5*peak2peak(double(pull(Scut)))';
            refsTmp = pull(Scut,'ref');
            lcut = length(ampsTmp);
            if strcmp(Scut(j).kstnm,"REVN")
                revn_amps(si:si+lcut-1) = ampsTmp;
                revn_t(si:si+lcut-1) = refsTmp;
            else
                revs_amps(si:si+lcut-1) = ampsTmp;
                revs_t(si:si+lcut-1) = refsTmp;
            end
        catch
            fprintf('couldnt strip rev of its amplitudes, continuing\n');
            continue;
        end
    end
    si = ei + 1;
end




