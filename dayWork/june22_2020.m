kstnms = ["BRUN";"BBIL";"BMAS";"BULB";"BPAT";"BRTU";"COHC";"CHSH";"PORT";"JSCH";"BOSC";"PKYU";"BREF";"PIAT";"PIS1"];
chans = ["BHZ";"BHZ";"BHZ";"BHZ";"BHZ";"HHZ";"HHZ";"HHZ";"HHZ";"HHZ";"HHZ";"HHZ";"BHZ";"HHZ";"HHZ"];

%%
maxN = 1200;
portTemplates = NaN(maxN,NtestTemplates);
jschTemplates = portTemplates;
boscTemplates = portTemplates;
pkyuTemplates = portTemplates;
brefTemplates = portTemplates;
piatTemplates = portTemplates;

%%
tw = 0.1;
noiseWin = 5;
secDur = 150;
newFs = 8;
lfc = 0.6;
hfc = 1.2;

%%
for i = 9:14 % length(kstnms)
    disp(' ');
    kstnm_ = kstnms(i);
    chan_ = chans(i);
    disp(kstnm_);
    disp(' ');
    
    %%
    Sorig = extractWaveforms(t-seconds(noiseWin),seconds(secDur),kstnm_,chan_);
    
    %%
    for j = 1:NtestTemplates
        disp(j)
        famI_ = newFamilies2{j};
        rawShifts_ = rawShifts{j};
        
        %%
        S = Sorig(famI_);
        S = differentiateWaveforms(S);
        S = detrendWaveforms(S);
        S = taperWaveforms(S,tw);
        S = filterWaveforms(S,lfc,hfc);
        S = taperWaveforms(S,tw);
        S = intWaveforms(S);
        S = detrendWaveforms(S);
        S = taperWaveforms(S,tw);
        S = resampleWaveforms(S,newFs);
        
        %%
        S = double(pull(S));
        S = apply_shifts(S,rawShifts_);
        S = normalizeWaveforms(S);
        S = pws(S);
        S = normalizeWaveforms(S);
        
        %%
        if i == 9
            portTemplates(:,j) = S(1:maxN);
        elseif i == 10
            jschTemplates(:,j) = S(1:maxN);
        elseif i == 11
            boscTemplates(:,j) = S(1:maxN);
        elseif i == 12
            pkyuTemplates(:,j) = S(1:maxN);
        elseif i == 13
            brefTemplates(:,j) = S(1:maxN);
        else
            piatTemplates(:,j) = S(1:maxN);
        end
    end
    
    %%
    clear Sorig;
end