function techniques = rsamTechniqueStr(meanFlags,rmsFlags)
ltechniques = length(meanFlags);
techniques = repmat("",ltechniques,1);
for i = 1:ltechniques
    meanFlag_ = meanFlags(i);
    rmsFlag_ = rmsFlags(i);
    if meanFlag_
        if rmsFlag_
            technique = "RMnS";
        else
            technique = "MnAm";
        end
    else
        if rmsFlag_
            technique = "RMdS";
        else
            technique = "MdAm";
        end
    end
    techniques(i) = technique;
end