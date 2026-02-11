clear; clear all; close all; clc;
%cd /home/shernandez/products/duplicates
%load ispt2016locatedtemplates.mat
cd ~/research/now/plata_repeaters/
load isptDataForPicking4.mat
lfc = 0.75; hfc = 6;
% Zf = filterWaveforms(differentiateWaveforms(table2struct(Z)),lfc,hfc);
% Nf = filterWaveforms(differentiateWaveforms(table2struct(N)),lfc,hfc);
% Ef = filterWaveforms(differentiateWaveforms(table2struct(E)),lfc,hfc);
Zf = table2struct(Z);
Nf = table2struct(N);
Ef = table2struct(E);

Zf = resampleWaveforms(Zf,40);
Nf = resampleWaveforms(Nf,40);
Ef = resampleWaveforms(Ef,40);
lT = length(Zf);
kN = 1;

warning off signal:findpeaks:largeMinPeakHeight

%%
pre = 0;
post = 15;
days = datetime(2017,01,144):datetime(2017,12,31);
lDays = length(days);
pks = cell(lDays,lT);
locs = pks;
mads = pks;

for i = 1:lDays
    tStart = days(i); %datetime(2016,01,i);
    disp(tStart);
    [yyyy_,mm_,ddd_] = datevec(tStart);
    if yyyy_ < 2017
        S(1) = loadWaveforms(tStart,1,Z.kstnm(1),Z.kcmpnm(1),Z.knetwk(1),Z.khole(1));
        %S(2) = loadWaveforms(tStart,1,N.kstnm(1),N.kcmpnm(1),N.knetwk(1),N.khole(1));
        %S(3) = loadWaveforms(tStart,1,E.kstnm(1),E.kcmpnm(1),E.knetwk(1),E.khole(1));
    else
        S(1) = loadWaveforms(tStart,1,Z.kstnm(1),Z.kcmpnm(1),Z.knetwk(1),"");
        %S(2) = loadWaveforms(tStart,1,N.kstnm(1),N.kcmpnm(1),N.knetwk(1),"");
        %S(3) = loadWaveforms(tStart,1,E.kstnm(1),E.kcmpnm(1),E.knetwk(1),"");
    end
    
    if any(isnat(pull(S,'ref')))
        continue
    else
        S = intWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc));
        S = resampleWaveforms(S,40);
        S = syncWaveforms(S);
        data = double(pull(S));
        winlen = 1 + ((post - pre)/S(1).delta);
        data2 = double(data).^2;
        box = ones(winlen,length(S));
        norms = fftfilt(box,data2);
        norms = sqrt(abs(norms));
        t = getTimeVec(S(1));
        
        parfor j = 1:lT
            M = Zf(j); %[Zf(j); Nf(j); Ef(j)];
            for k = 1:kN%:3
                M(k) = cutWaveforms(M(k),M(k).ref+seconds(3),pre,post);
            end
            M = syncWaveforms(M);
            M = normalizeWaveforms(M,true,true);
            
            templates = normalizeWaveforms(pull(M));
            winlen = size(templates,1);
            ccnorm = fftfilt(double(flipud(templates)),data);
            snorms = norms > 1;
            ccnorm = snorms.*ccnorm./norms;
            
            stack = nanmean(ccnorm,2); mad_ = nanstd(stack);
            newThresh_ = min([0.7 15*mad_]);
            [pks_,locs_] = findpeaks(stack,'MinPeakHeight',newThresh_,'MinPeakDistance',winlen);
            pks{i,j} = pks_;
            mads{i,j} = pks_/mad_;
            locs{i,j} = datenum(t(locs_));
        end
    end
end

%%
tI = [];
tMaster = [];
pksMaster = [];
for i = 1:lT
    tnew = dn2dt(cat(1,locs{:,i}));
    tMaster = [tMaster; tnew];
    pksMaster = [pksMaster; cat(1,pks{:,i})];
    tI = [tI; i*ones(length(tnew),1)];
end
[tMaster,sI] = sort(tMaster);
pksMaster = pksMaster(sI);
tI = tI(sI);


[t_,cc_,removeIndices] = removeRepeatedMatches(tMaster,pksMaster,30);
for jj = 1:length(removeIndices)
    tI(removeIndices{jj}) = [];
end

ccI = cc_ >= 0.5;
t_ = t_(ccI);
cc_ = cc_(ccI);
tI = tI(ccI);

close all
figure(); S = scatter(t_,1:length(t_),[],cc_,'o','filled'); zoom on; colorbar;