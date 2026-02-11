function cc = stretch(refTrace,testTrace,base,searchVec,si,ei,detrendFlag)
lstretches = length(searchVec);
tc = size(testTrace,2);
cc = NaN(lstretches,tc);

for i = 1:lstretches
    refTrace_ = resample(refTrace,searchVec(i),base);
    refTrace_ = refTrace_(si:ei,:);
    refTrace_ = normalizeWaveforms(refTrace_,detrendFlag,0);
    cc_ = sum(refTrace_.*testTrace);
    cc(i,:) = cc_;
end