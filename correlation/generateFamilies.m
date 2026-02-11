function [multipletI,lPerMultiplet,nSingletons,multipletNumber,tree] = generateFamilies(maxccp,...
    thresh,method,minNumberOfMembers,plotFlag)
if nargin < 5; plotFlag = false; end
if nargin < 4; minNumberOfMembers = 2; end
if nargin < 3; thresh = 0.8; end
if nargin < 2; method = "complete"; end

%%
if ~isvector(maxccp)
    fprintf(2,"Error: D is not a vector. Doing nothing.\n");
    return;
end

%%
[multipletI,lPerMultiplet,multipletNumber,tree] = returnClusts(1-maxccp',thresh,method);
goodI = lPerMultiplet >= minNumberOfMembers;
nClusteredEvents = sum(lPerMultiplet(goodI));
nMultiplets = sum(goodI);
nSingletons = sum(lPerMultiplet(~goodI)); %number unclustered

%%
if ~plotFlag
    return;
end

%%
mostMembers = lPerMultiplet(1);
markersize = 20;
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
plot(1:nMultiplets,lPerMultiplet(1:nMultiplets),'.','markersize',markersize);
xlabel('multiplet number');
titleStr = sprintf("number of members per multiplet (min. of %d members); Method: %s",...
    minNumberOfMembers,method);
title(titleStr);
text(0.8*nMultiplets,0.95*mostMembers,...
    sprintf("number of singletons: %d",nSingletons));
text(0.78*nMultiplets,0.93*mostMembers,...
    sprintf("number of clustered events: %d",nClusteredEvents));
zoom on;