function [newEvents,newFamilies,lPerMultiplet,...
    singletonI,lSingletons,maxccp] = ...
    pruneAndMergeEvents(indiv_events,...
    maxccpOrig,threshVec,NMINVEC,method,PWSFLAG,CANUSEGPU,MAXLAG,NCHAN)
%
% extensive rewrites 08 DECEMBER 2025
%

%%
if nargin < 9; NCHAN = 1; end
if nargin < 8; MAXLAG = []; end
if nargin < 7; CANUSEGPU = false; end
if nargin < 6; PWSFLAG = false; end
if nargin < 5; method = "complete"; end
if nargin < 4; NMINVEC = 8; end
if nargin < 3; threshVec = 0.8; end

if sum(~isfinite(indiv_events))
    fprintf("there are events that have non-finites in them\n");
    return;
end

%%
lThreshes = length(threshVec);
thisThresh = threshVec(1);
NMIN = NMINVEC(1);

%%
indiv_events = normalizeWaveforms(indiv_events);
[winlen,nEvents] = size(indiv_events);
if isempty(MAXLAG)
    MAXLAG = winlen - 1;
end
disp(table(method,PWSFLAG,CANUSEGPU,MAXLAG,NCHAN))
fprintf("\n");

%%
fprintf("iteration: 1; nmin: %d; thresh: %f\n",NMIN,thisThresh);
[multipletI,lPerMultiplet] = generateFamilies(maxccpOrig,thisThresh,...
    method,2,false);
maxccpOrig = squareform(maxccpOrig);
goodI = lPerMultiplet >= NMIN;
nClusteredEvents = sum(lPerMultiplet(goodI));
nUnclustered = sum(lPerMultiplet(~goodI));
nStacks = sum(goodI);
nEventsNew = nUnclustered + nStacks;
lPerMultiplet = lPerMultiplet(goodI);

if ~nStacks
    fprintf("no families found for threshold level: %f\n",thisThresh);
    return;
end

singletonI = true(nEvents,1);
clusteredI = cat(1,multipletI{1:nStacks});
singletonI(clusteredI) = false;
lSingletons = sum(singletonI); %nUnclustered; %length(singletonIndexes);
singletonI = find(singletonI);

newFamilies = multipletI(goodI);
newEvents = indiv_events(:,1:nEventsNew); %simple preallocation
newEvents(:,nStacks+1:end) = indiv_events(:,singletonI); %populate the singletons (the stacks will come later)

%% be verbose
fprintf("%d event(s) being compressed into %d stack(s)\n",nClusteredEvents,nStacks);
fprintf("number of isolated/unpaired events: %d\n",lSingletons);
fprintf("new reduced number of events: %d\n",nEventsNew);
fprintf("number of events pruned as probably redundant: %d\n",nClusteredEvents - nStacks);

fprintf("generating stacks for each multiplet\n");
chan_length = winlen/NCHAN;
if chan_length ~= round(chan_length)
    fprintf("wrong number of channels specified. Double check the actual value\n");
    return;
end
for j = 1:nStacks
    fam_ = newFamilies{j}; %fam_ is a regular (normal) array, no funny stuff (not a cell)
    indivTmp = normalizeWaveforms(detrend(indiv_events(:,fam_)));
    indivTmp = apply_vdcc(indivTmp,[],false,false,false,CANUSEGPU);% align and stack waveforms to generate new stack
    lFam = length(fam_);
    indivTmp = reshape(indivTmp,[],NCHAN,lFam);
    indivTmp = detrend(indivTmp);
    indivTmp = indivTmp./rssq(indivTmp);

    thisStack = zeros(chan_length,NCHAN);
    for k = 1:NCHAN
        indivTmp_ = squeeze(indivTmp(:,k,:));
        if PWSFLAG
            thisStack(:,k) = normalizeWaveforms(pws(indivTmp_));
        else
            [U,~] = svd(indivTmp_,"econ");
            U = U(:,1);
            mean_pol = sign(mean(sign(sum(U.*indivTmp_))));
            thisStack(:,k) = normalizeWaveforms(mean_pol*U);
            %newEvents(:,j) = normalizeWaveforms(mean(normalizeWaveforms(detrend(indivTmp)),2,"omitnan")); % stick the shifted/aligned stacks where they belong
        end
    end
    if NCHAN > 1
        thisStack = normalizeWaveforms(thisStack(:));
    end
    newEvents(:,j) = thisStack;
end

%% here we do new cross correlations using vectors from new stacks
fprintf("cross correlating new stacks with older singletons\n");
timer_val = tic;
refBlock = normalizeWaveforms(newEvents(:,1:nStacks));
dSingletons = normalizeWaveforms(indiv_events(:,singletonI));
if CANUSEGPU
    refBlock = gpuArray(single(refBlock));
    dSingletons = gpuArray(single(dSingletons));
end
maxabsccp = zeros(nEventsNew);
[maxccpUneven,~] = ...
    doccFreqCircShiftPolarities(dSingletons,false,refBlock,MAXLAG);
maxabsccp(1:nStacks,:) = maxccpUneven;
maxabsccp(nStacks+1:end,1:nStacks) = maxccpUneven(:,nStacks+1:end)';
maxabsccp(nStacks+1:end,nStacks+1:end) = maxccpOrig(singletonI,singletonI);
maxccp = squareform(maxabsccp)'; %always convert to column (vertical) vector please...
elapsed_time = toc(timer_val);
fprintf("elapsed time iteration 1: %f\n",elapsed_time);

nMultiplets = length(newFamilies);
newFamilies = [newFamilies(:); cell(lSingletons,1)];
for i = 1:lSingletons
    newFamilies{nMultiplets+i} = singletonI(i);
end

if lThreshes < 2
    return;
end

nStacksOrig = nStacks;
singletonIorig = singletonI;
for i = 2:lThreshes
    fprintf("\n");
    thisThresh = threshVec(i);
    NMIN = NMINVEC(i);
    fprintf("iteration: %d/%d; nmin: %d; thresh: %f\n",i,lThreshes,NMIN,thisThresh);
    [multipletI,lPerMultiplet] = generateFamilies(maxccp,thisThresh,...
        method,2,false);
    goodI = lPerMultiplet >= 2;
    nStacks = sum(goodI);
    if ~nStacks
        fprintf("no families found for threshold level: %f\n",thisThresh);
        return;
    end
    delI = false(nStacksOrig,1);
    new_multiplet_index = nStacksOrig;
    for j = 1:nStacks
        fam_ = multipletI{j};
        mergeI = fam_ <= nStacksOrig;
        dupsExist = sum(mergeI);
        new_members = singletonIorig(fam_(~mergeI)-nStacksOrig);
        if dupsExist
            dupI = fam_(mergeI);
            dupI_ = dupI(1);
            newFamilies{dupI_} = unique([cat(1,newFamilies{dupI});...
                new_members]);
            if dupsExist > 1
                delI(dupI(2:end)) = true; %answer to question: yes, delete
            end
        else
            new_multiplet_index = new_multiplet_index + 1;
            newFamilies{new_multiplet_index} = new_members;
        end
    end
    newFamilies = newFamilies(1:new_multiplet_index);
    newFamilies(delI) = []; %answer to question: yes, delete
    nStacks = length(newFamilies);
    lPerMultiplet = ones(nStacks,1);
    for j = 1:nStacks
        famI = newFamilies{j};
        lPerMultiplet(j) = length(famI);
    end
    [lPerMultiplet,sI] = sort(lPerMultiplet,"descend");
    newFamilies = newFamilies(sI);
    singletonI = true(nEvents,1);
    clusteredI = cat(1,newFamilies{:});
    singletonI(clusteredI) = false;

    goodI = lPerMultiplet >= NMIN;
    nStacks = sum(goodI);
    if ~nStacks
        fprintf("no families found for threshold level: %f\n",thisThresh);
        return;
    end

    nClusteredEvents = sum(lPerMultiplet(goodI));
    lSingletons = sum(singletonI);
    nUnclustered = sum(lPerMultiplet(~goodI)) + lSingletons;
    nEventsNew = nUnclustered + nStacks;
    lPerMultiplet = lPerMultiplet(goodI);
    clusteredI = cat(1,newFamilies{~goodI}); %stranded....
    singletonI(clusteredI) = true;
    lSingletons = sum(singletonI); %nUnclustered; %length(singletonIndexes);
    singletonI = find(singletonI); %returns sorted

    newFamilies = newFamilies(goodI);
    newEvents = indiv_events(:,1:nEventsNew); %simple preallocation
    newEvents(:,nStacks+1:end) = indiv_events(:,singletonI);

    %% be verbose
    fprintf("%d event(s) being compressed into %d stack(s)\n",nClusteredEvents,nStacks);
    fprintf("number of isolated/unpaired events: %d\n",lSingletons);
    fprintf("new reduced number of events: %d\n",nEventsNew);
    fprintf("number of events pruned as probably redundant: %d\n",nClusteredEvents - nStacks);

    %%
    fprintf("generating stacks for each multiplet\n");
    for j = 1:nStacks
        fam_ = newFamilies{j}; %fam_ is a regular (normal) array, no funny stuff (not a cell)
        indivTmp = normalizeWaveforms(detrend(indiv_events(:,fam_)));
        indivTmp = apply_vdcc(indivTmp,[],false,false,false,CANUSEGPU);% align and stack waveforms to generate new stack
        lFam = length(fam_);
        indivTmp = reshape(indivTmp,[],NCHAN,lFam);
        indivTmp = detrend(indivTmp);
        indivTmp = indivTmp./rssq(indivTmp);

        thisStack = zeros(chan_length,NCHAN);
        for k = 1:NCHAN
            indivTmp_ = squeeze(indivTmp(:,k,:));
            if PWSFLAG
                thisStack(:,k) = normalizeWaveforms(pws(indivTmp_));
            else
                [U,~] = svd(indivTmp_,"econ");
                U = U(:,1);
                mean_pol = sign(mean(sign(sum(U.*indivTmp_))));
                thisStack(:,k) = normalizeWaveforms(mean_pol*U);
                %newEvents(:,j) = normalizeWaveforms(mean(normalizeWaveforms(detrend(indivTmp)),2,"omitnan")); % stick the shifted/aligned stacks where they belong
            end
        end
        if NCHAN > 1
            thisStack = normalizeWaveforms(thisStack(:));
        end
        newEvents(:,j) = thisStack;
    end

    fprintf("cross correlating new stacks with older singletons\n");
    timer_val = tic;
    refBlock = normalizeWaveforms(newEvents(:,1:nStacks));
    dSingletons = normalizeWaveforms(indiv_events(:,singletonI));
    if CANUSEGPU
        refBlock = gpuArray(single(refBlock));
        dSingletons = gpuArray(single(dSingletons));
    end
    maxabsccp = zeros(nEventsNew);
    [maxccpUneven,~] = ...
        doccFreqCircShiftPolarities(dSingletons,false,refBlock,MAXLAG);
    maxabsccp(1:nStacks,:) = maxccpUneven;
    maxabsccp(nStacks+1:end,1:nStacks) = maxccpUneven(:,nStacks+1:end)';
    maxabsccp(nStacks+1:end,nStacks+1:end) = maxccpOrig(singletonI,singletonI);
    maxccp = squareform(maxabsccp)'; %always convert to column (vertical) vector please...
    elapsed_time = toc(timer_val);
    fprintf("elapsed time iteration %d: %f\n",i,elapsed_time);

    %%
    nStacksOrig = nStacks;
    singletonIorig = singletonI;
end

%%
nMultiplets = length(newFamilies);
newFamilies = [newFamilies(:); cell(lSingletons,1)];
for i = 1:lSingletons
    newFamilies{nMultiplets+i} = singletonI(i);
end