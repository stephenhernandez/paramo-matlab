function [travTimes,dists,depths] = getTravelTimes(velModel,phaseList,depths,dists,vectorFlag)

if nargin < 1; velModel = []; end
if nargin < 2; phaseList = 'p,P,Pdiff'; end
if nargin < 3; depths = 0.5 + (0:799)'; end
if nargin < 4; dists = 0.5 + (0:10:20000)'; end
if nargin < 5; vectorFlag = false; end

%%
lDepths = length(depths);
lDists = length(dists);

if vectorFlag
    if lDepths ~= lDists
        fprintf('length of dists and depth vectors are not the same, exiting\n');
        return;
    end
    travTimes = NaN(lDepths,1);
    for i = 1:lDepths
        depth = depths(i);
        dist = dists(i);
        tt = taupTime(velModel,depth,phaseList,'km',dist);
        if ~isempty(tt)
            travTimes(i) = tt(1).time;
        end
    end
else
    travTimes = NaN(lDepths,lDists);
    for i = 1:lDepths
        disp(i);
        depth = depths(i);
        parfor j = 1:lDists
            dist = dists(j);
            tt = taupTime(velModel,depth,phaseList,'km',dist);
            if ~isempty(tt)
                travTimes(i,j) = tt(1).time;
            end
        end
    end
end