function [dVcaus,dVacaus,dVsym,ccCaus,ccAcaus,ccSym,...
    rmsPrcntCaus,rmsPrcntAcaus,rmsPrcntSym] = ...
    pii(t,caus,acaus,symStack,referenceStartTime,cwiStart,cwiEnd,newFs,...
    stnm1,stnm2,chan1,chan2,maxPrcnt,stackMethod,technique,force32bit,lfc,hfc)

%
% passive image interferometry
%


%%
plotFlag = false;

%%
if nargin < 14
    stackMethod = 'mean';
end

if nargin < 15
    technique = 'stretch';
end

if nargin < 16
    force32bit = false;
end

%%
stretchFlag = false;
if strcmp(technique,'stretch')
    stretchFlag = true;
end

%%
if strcmp(stackMethod,'pws')
    refTrace = pws(caus(:,t <= referenceStartTime));
elseif strcmp(stackMethod,'med')
    refTrace = nanmedian(caus(:,t <= referenceStartTime),2);
elseif strcmp(stackMethod,'mean')
    refTrace = nanmean(caus(:,t <= referenceStartTime),2);
end
testTrace = caus;

if stretchFlag
    [dVcaus,ccCaus,~,~,rmsPrcntCaus] = cwiShifts(refTrace,testTrace,cwiStart,cwiEnd,newFs,plotFlag,maxPrcnt,lfc,hfc);
end

%%
if ~strcmp(chan1,chan2) || ~strcmp(stnm1,stnm2)
    if strcmp(stackMethod,'pws')
        refTrace = pws(acaus(:,t <= referenceStartTime));
    elseif strcmp(stackMethod,'med')
        refTrace = nanmedian(acaus(:,t <= referenceStartTime),2);
    elseif strcmp(stackMethod,'mean')
        refTrace = nanmean(acaus(:,t <= referenceStartTime),2);
    end
    testTrace = acaus;
    
    
    if stretchFlag
        [dVacaus,ccAcaus,~,~,rmsPrcntAcaus] = cwiShifts(refTrace,testTrace,cwiStart,cwiEnd,newFs,plotFlag,maxPrcnt,lfc,hfc);
    end
else
    dVacaus = dVcaus;
    ccAcaus = ccCaus;
    rmsPrcntAcaus = rmsPrcntCaus;
end

%%
if ~strcmp(chan1,chan2) || ~strcmp(stnm1,stnm2)
    if strcmp(stackMethod,'pws')
        refTrace = pws(symStack(:,t <= referenceStartTime));
    elseif strcmp(stackMethod,'med')
        refTrace = nanmedian(symStack(:,t <= referenceStartTime),2);
    elseif strcmp(stackMethod,'mean')
        refTrace = nanmean(symStack(:,t <= referenceStartTime),2);
    end
    
    testTrace = symStack;
    
    if stretchFlag
        [dVsym,ccSym,~,~,rmsPrcntSym] = cwiShifts(refTrace,testTrace,cwiStart,cwiEnd,newFs,plotFlag,maxPrcnt,lfc,hfc);
    end
else
    dVsym = dVcaus;
    ccSym = ccCaus;
    rmsPrcntSym = rmsPrcntCaus;
end

%%
if force32bit
    dVcaus = single(dVcaus);
    dVacaus = single(dVacaus);
    dVsym = single(dVsym);
    ccCaus = single(ccCaus);
    ccAcaus = single(ccAcaus);
    ccSym = single(ccSym);
end