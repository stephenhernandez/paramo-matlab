function [Zeros,Poles,Constant,Gain,Sensitivity,Normalization,I,ic,...
    Stla,Stlo,Stel,CMPNM,instrumentType] = getPolesZeros(S,t)
if nargin < 2
    t = pull(S,"ref");
end

%%
lT = length(t);
load('~/igdata/ecuadorSensorDataTable10','kstnm','knetwk','kcmpnm','khole',...
    'ZPT','tStart','constant','sensitivity','a0','stla','stlo','stel');

instType = ZPT.instType;
instGain = ZPT.instGain;
kcmpnms = pull(S,"kcmpnm");
kholes = pull(S,"khole");
kstnms = pull(S,"kstnm");
knetwks = pull(S,"knetwk");
sncl = strcat(knetwks,kstnms,kholes,kcmpnms);

%%
[uniqueSNCLs,ia,ic] = unique(sncl);
luniqueSNCLs = length(uniqueSNCLs);
uniqkcmpnms = kcmpnms(ia);
uniqkholes = kholes(ia);
uniqkstnms = kstnms(ia);
uniqknetwks = knetwks(ia);

Zeros = cell(lT,1);
Poles = Zeros;
Constant = NaN(lT,1);
Gain = Constant;
Sensitivity = Constant;
Normalization = Constant;
I = Constant;
Stla = Constant;
Stlo = Stla;
Stel = Stla;
CMPNM = repmat("",lT,1);
instrumentType = CMPNM;

%%
for j = 1:luniqueSNCLs
    ui = find(j == ic);
    kstnm_ = uniqkstnms(j);
    knetwk_ = uniqknetwks(j);
    kcmpnm_ = uniqkcmpnms(j);
    khole_ = uniqkholes(j);

    %%
    lia = ismember(kstnm,kstnm_) & ismember(knetwk,knetwk_) & ...
        ismember(kcmpnm,kcmpnm_) & ismember(khole,khole_);
    if ~sum(lia)
        fprintf("something went wrong, cannot find index to SNCL: %s\n",uniqueSNCLs(j));
        return;
    end
    liaOrig = find(lia);
    lSNCLRepeats = length(ui);
    tStart_ = tStart(liaOrig);
    for i = 1:lSNCLRepeats
        t_ = t(ui(i));
        tI = find(t_ >= tStart_);
        if isempty(tI)
            tI = 1;
        end
        tI = tI(end);

        %%
        lia = liaOrig(tI);
        if lT > 1
            Zeros{ui(i)} = ZPT.zeros{lia};
            Poles{ui(i)} = ZPT.poles{lia};
        else
            Zeros = ZPT.zeros{lia};
            Poles = ZPT.poles{lia};
        end

        %%
        Constant(ui(i)) = constant(lia);
        Sensitivity(ui(i)) = sensitivity(lia);
        Normalization(ui(i)) = a0(lia);
        I(ui(i)) = tI;
        Stla(ui(i)) = stla(lia);
        Stlo(ui(i)) = stlo(lia);
        Stel(ui(i)) = stel(lia);

        %%
        gain_ = instGain(lia);
        if iscell(gain_)
            gain_ = gain_{1};
        end
        Gain(ui(i)) = gain_;
        CMPNM(ui(i)) = kcmpnm_;
        instrumentType(ui(i)) = instType(lia);
    end
end