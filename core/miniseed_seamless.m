function trace = miniseed_seamless(trace,bmI,cBlockMatrix,blockMatrix,x0,...
    numelColumn,si)

for i = 2:cBlockMatrix
    numelColumn_ = numelColumn(i);
    bmI_ = bmI(:,i);
    column = blockMatrix(bmI_,i);
    x0_ = x0(i);
    dNext = cumsum([x0_; column(2:numelColumn_)]);
    trace(si:si+numelColumn_-1,1) = dNext;
    si = si + numelColumn_;
end

