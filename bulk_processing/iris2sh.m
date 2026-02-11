function S = iris2sh(mytrace)
%%
lm = length(mytrace);
if lm < 1
    S = populateWaveforms(1);
    fprintf('no iris data to convert\n');
    return;
end

%%
S = populateWaveforms(lm);

%%
for i = 1:lm
    knetwk = string(mytrace(i).network);
    kstnm = string(mytrace(i).station);
    khole = string(mytrace(i).location);
    kcmpnm = string(mytrace(i).channel);
    stla = mytrace(i).latitude;
    stlo = mytrace(i).longitude;
    d = mytrace(i).data;
    npts = mytrace(i).sampleCount;
    delta = 1./mytrace(i).sampleRate;
    ref = mytrace(i).startTime;
    az = mytrace(i).azimuth;
    inc = mytrace(i).dip;
    
    S(i).ref = dn2dt(ref);
    S(i).b = 0;
    S(i).kstnm = kstnm;
    S(i).khole = khole;
    S(i).kcmpnm = kcmpnm;
    S(i).knetwk = knetwk;
    S(i).stla = stla;
    S(i).stlo = stlo;
    S(i).npts = npts;
    S(i).delta = delta;
    S(i).d = d;
    S(i).cmpinc = inc;
    S(i).cmpaz = az;
    S(i).e = seconds((S(i).npts - 1)*S(i).delta);
end