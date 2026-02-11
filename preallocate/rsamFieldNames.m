function fieldName = rsamFieldNames(lfc,hfc,technique)

%%
if lfc < 1
    secs = 1/lfc;
    if isflint(secs)
        f1 = sprintf("%dSec",secs);
    else
        f1 = sprintf("%3.2fSec",secs);
    end
else
    if isflint(lfc)
        f1 = sprintf("%dHz",lfc);
    else
        f1 = sprintf("%3.2Hz",lfc);
    end
end

%%
if hfc < 1
    secs = 1/hfc;
    if isflint(secs)
        f2 = sprintf("%dSec",secs);
    else
        f2 = sprintf("%3.2fSec",secs);
    end
else
    if isflint(hfc)
        f2 = sprintf("%dHz",hfc);
    else
        f2 = sprintf("%3.2fHz",hfc);
    end
end

if ~isfinite(hfc)
    f2 = "";
end

%%
fieldName = strcat(technique,"_",f1,f2);
fieldName = matlab.lang.makeValidName(fieldName);
