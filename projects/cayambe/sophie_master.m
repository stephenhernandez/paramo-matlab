%cd ~/research/now/cayambe/pedernales_triggering/
t = [datetime(2016,09,05,13,58,59.010000);...
    datetime(2016,12,13,05,16,26.009994);...
    datetime(2016,06,09,22,03,09.200006)];

%%
newFs = 100;
noiseWinLen = 15;
noiseNpts = noiseWinLen*newFs;
sigLen = noiseNpts;

nfft = 2^nextpow2(noiseNpts);

%%
tw = 0.02;
close all;
for i = 1:length(t)
    legStr = [];
    S(1) = extractWaveforms(t(i)-seconds(20),seconds(40),"ANGU","SHZ","EC");
    S(2) = extractWaveforms(t(i)-seconds(20),seconds(40),"CAYA","HHZ","EC");
    S(3) = extractWaveforms(t(i)-seconds(20),seconds(40),"CAYR","SHZ","EC");
    S = detrendWaveforms(S);
    S = interpolateWaveforms(S);
    S = resampleWaveforms(S,newFs);
    %plotWaveforms(S);
    refs = pull(S,'ref');
    sumgood = sum(~isnat(refs));
    S_ = S(~isnat(refs));
    for j = 1:sumgood
        S__ = S_(j);
        d = S__.d;
        dnoise = taper(detrend(d(1:noiseNpts)),tw);
        dsig = taper(detrend(d(noiseNpts+1:noiseNpts+1+noiseNpts-1)),tw);
        [pxxN_,fxx] = pmtm(dnoise,{2,'trace'},nfft,newFs);
        pxxS_ = pmtm(dsig,{2,'trace'},nfft,newFs);
        figure(i); hold on;
        pp = loglog(fxx(2:end),pxxS_(2:end)./pxxN_(2:end),'-','linewidth',3); zoom on; grid on;
        pp.Color(4) = 0.8;
        legStr = [legStr; string(S__.kstnm)];
        %title([num2str(i),': ',S__.kstnm]);
    end
    legend(legStr);
end

%%
for i = 1:length(t)
    legStr = [];
    S(1) = extractWaveforms(t(i)-seconds(20),seconds(40),"ANGU","SHZ","EC");
    S(2) = extractWaveforms(t(i)-seconds(20),seconds(40),"CAYA","HHZ","EC");
    S(3) = extractWaveforms(t(i)-seconds(20),seconds(40),"CAYR","SHZ","EC");
    S = detrendWaveforms(S);
    S = interpolateWaveforms(S);
    S = resampleWaveforms(S,newFs);
    plotWaveforms(S,6,16);
end