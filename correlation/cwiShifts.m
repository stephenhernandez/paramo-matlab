function [prcntChange,peakCC,cc,err_rms_prcnt] = cwiShifts(refTrace,testTrace,...
    secStart,secEnd,inputFs,plotFlag,maxPrcnt,lfc,hfc,detrendFlag)
if nargin < 3; secStart = 3; end
if nargin < 4; secEnd = 12; end
if nargin < 5; inputFs = 100; end
if nargin < 6; plotFlag = false; end
if nargin < 7; maxPrcnt = 3; end
if nargin < 10; detrendFlag = true; end

%%
si = 1+floor(secStart*inputFs);
ei = 0+floor(secEnd*inputFs);

%%
testTrace = testTrace(si:ei,:);
testTrace = normalizeWaveforms(testTrace,detrendFlag);
rc = size(refTrace,2);
if rc > 1
    fprintf("will run a little slower since there is a unique ref trace for each test trace\n");
end

%% granularity: 1%
fprintf("granularity: 1%%\n");
base = 100;
prcntInc = round(maxPrcnt*base/100);
epsSearch1 = (base-prcntInc:base+prcntInc)'; %inc. 1%, but only +/- maxPrcnt
cc = stretch(refTrace,testTrace,base,epsSearch1,si,ei,detrendFlag);
ccI = cc <= 0;
cc(ccI) = 0;
[~,maxccI] = max(cc,[],1); % argmaxx
sampleShiftsThisBase = epsSearch1(maxccI);
[~,maxccIlinear] = max(cc,[],1,"linear"); % argminx
peakCC = cc(maxccIlinear);

%% granularity: 0.1%
fprintf("granularity: 0.1%%\n");
prevBase = base;    % old base
base = 1000;        % new base
epsSearch01 = (base-10:base+10)'; %this is relative to the peak so should never change
uniqueShifts = unique(sampleShiftsThisBase)';
lI = length(uniqueShifts);

for i = 1:lI
    thisShift = uniqueShifts(i);
    thisBlock = sampleShiftsThisBase == thisShift;
    %fprintf("%d %d %d\n",i,sum(thisBlock),lI);
    thisTestTrace = testTrace(:,thisBlock);
    if rc > 1
        refTrace_ = resample(refTrace(:,thisBlock),thisShift,prevBase);
    else
        refTrace_ = resample(refTrace(:,1),thisShift,prevBase);
    end

    cc_ = stretch(refTrace_,thisTestTrace,base,epsSearch01,si,ei,detrendFlag); %ref trace is close to full length (havent cut yet)
    ccI = cc_ <= 0;
    cc_(ccI) = 0;
    [~,maxccI] = max(cc_,[],1);
    delta_prcntShifts = epsSearch01(maxccI);
    sampleShiftsThisBase(thisBlock) = (delta_prcntShifts-base) + base*thisShift/prevBase;
    [~,maxccIlinear] = max(cc_,[],1,"linear");
    peakCC(thisBlock) = cc_(maxccIlinear);
end

%% granularity: 0.01%
fprintf("granularity: 0.01%%\n");
prevBase = base;    % old base
base = 10000;        % new base
epsSearch001 = (base-10:base+10)'; %this is relative to the peak so should never change
uniqueShifts = unique(sampleShiftsThisBase)';
lI = length(uniqueShifts);

for i = 1:lI
    thisShift = uniqueShifts(i);
    thisBlock = sampleShiftsThisBase == thisShift;
    %fprintf("%d %d %d\n",i,sum(thisBlock),lI);
    thisTestTrace = testTrace(:,thisBlock);
    if rc > 1
        refTrace_ = resample(refTrace(:,thisBlock),thisShift,prevBase);
    else
        refTrace_ = resample(refTrace(:,1),thisShift,prevBase);
    end

    cc_ = stretch(refTrace_,thisTestTrace,base,epsSearch001,si,ei,detrendFlag); %ref trace is close to full length (havent cut yet)
    ccI = cc_ <= 0;
    cc_(ccI) = 0;
    [~,maxccI] = max(cc_,[],1);
    delta_prcntShifts = epsSearch001(maxccI);
    sampleShiftsThisBase(thisBlock) = (delta_prcntShifts-base) + base*thisShift/prevBase;
    [~,maxccIlinear] = max(cc_,[],1,"linear");
    peakCC(thisBlock) = cc_(maxccIlinear);
end

%% granularity: 0.0025%
fprintf("granularity: 0.0025%%\n");
prevBase = base;        % old base
base = 40000;           % new base
epsSearch00025 = (base-4:base+4)'; %this is relative to the peak so should never change
uniqueShifts = unique(sampleShiftsThisBase)';
lI = length(uniqueShifts);

for i = 1:lI
    thisShift = uniqueShifts(i);
    thisBlock = sampleShiftsThisBase == thisShift;
    %fprintf("%d %d %d\n",i,sum(thisBlock),lI);
    thisTestTrace = testTrace(:,thisBlock);
    if rc > 1
        refTrace_ = resample(refTrace(:,thisBlock),thisShift,prevBase);
    else
        refTrace_ = resample(refTrace(:,1),thisShift,prevBase);
    end

    cc_ = stretch(refTrace_,thisTestTrace,base,epsSearch00025,si,ei,detrendFlag); %ref trace is close to full length (havent cut yet)
    ccI = cc_ <= 0;
    cc_(ccI) = 0;
    [~,maxccI] = max(cc_,[],1);
    delta_prcntShifts = epsSearch00025(maxccI);
    sampleShiftsThisBase(thisBlock) = (delta_prcntShifts-base) + base*thisShift/prevBase;
    [~,maxccIlinear] = max(cc_,[],1,"linear");
    peakCC(thisBlock) = cc_(maxccIlinear);
end

%%
prcntChange = 100*(1-(sampleShiftsThisBase/base));
peakCC = peakCC';

if nargout > 3
    if ~exist('lfc','var') || ~exist('hfc','var')
        fprintf(2,'Error, lfc/hfc do not exist, cannot compute error bars, exiting...\n');
        return;
    end
    cf = sqrt(lfc*hfc);
    omega_c = 2*pi*cf;
    err_rms_prcnt = ...
        100*(sqrt((1-peakCC).*(1+peakCC))./...
        (2*peakCC)).*(sqrt((6*sqrt(pi/2)./...
        (hfc - lfc))./((omega_c.^2)*((secEnd^3) - (secStart^3)))));
end

%%
if ~plotFlag
    return;
end
figure(); plot(sampleShiftsThisBase,peakCC,'o'); zoom on;
figure(); plot(peakCC,'o'); zoom on;
figure(); plot(100*(base-eps)/base,cc,'.'); zoom on;