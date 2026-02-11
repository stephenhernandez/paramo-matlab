function S = duplicateStructDifferentDay(S,newDay)%,rawDataDir)
%if nargin < 5; rawDataDir = '~/rawdata'; end

lS = length(S);
for i = 1:lS
    kstnm_ = S(i).kstnm;
    kntwk_ = S(i).knetwk;
    kcmpnm_ = S(i).kcmpnm;
    khole_ = S(i).khole;
    S_ = loadWaveforms(newDay,1,kstnm_,kcmpnm_,kntwk_,khole_); %,true,true,rawDataDir);
    S(i) = S_;
end