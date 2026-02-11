function [tMain,ccMain,ampMain,dMag,magMain,templateNumber,...
    madMain,nUsedMain,idMain] = loadRepeaterCatalog2(dirName)
mainDir = fullfile('~','templateSearch',dirName);
cd(fullfile(mainDir));
files = dir('dayTemplateSearch*txt');
lfiles = length(files);

%%
nMax = 1e6;
tMain = NaT(nMax,1);
idMain = repmat("",nMax,1);
ccMain = NaN(nMax,1);
ampMain = ccMain;
dMag = ccMain;
magMain = ccMain;
templateNumber = ccMain;
madMain = ccMain;
nUsedMain = ccMain;

formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %s';
%%
n = 1;
for i = 1:lfiles
    fid = fopen(files(i).name);
    T = textscan(fid,formatSpec);
    t_ = datetime(cat(2,T{1:6}));
    lEvents = length(t_);
    fprintf("%s: %s, %d\n",mainDir,files(i).name,lEvents);

    cc_ = T{9}; %T.ccMain;
    amp_ = T{8}; %T.ampMain;
    dmag_ = T{7}; %T.dMag;
    tnum_ = T{10}; %T.templateNumber;
    mad_ = T{11}; %T.madMain;
    mag_ = T{12}; %T.magMain;
    nused_ = T{13}; %T.nUsedMain;
    id_ = T{14}; %T.templateNumber;

    tMain(n:n+lEvents-1) = t_;
    ccMain(n:n+lEvents-1) = cc_;
    idMain(n:n+lEvents-1) = id_;
    ampMain(n:n+lEvents-1) = amp_;
    dMag(n:n+lEvents-1) = dmag_;
    magMain(n:n+lEvents-1) = mag_;
    templateNumber(n:n+lEvents-1) = tnum_;
    madMain(n:n+lEvents-1) = mad_;
    nUsedMain(n:n+lEvents-1) = nused_;

    n = n + lEvents;
    fclose(fid);
end

%%
tMain(n:end) = [];
dMag(n:end) = [];
ampMain(n:end) = [];
ccMain(n:end) = [];
idMain(n:end) = [];
templateNumber(n:end) = [];
magMain(n:end) = [];
madMain(n:end) = [];
nUsedMain(n:end) = [];

%%
[tMain,sI] = sort(tMain);
dMag = dMag(sI);
ampMain = ampMain(sI);
ccMain = ccMain(sI);
templateNumber = templateNumber(sI);
magMain = magMain(sI);
madMain = madMain(sI);
nUsedMain = nUsedMain(sI);
idMain = idMain(sI);
