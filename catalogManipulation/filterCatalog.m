function [t,cc,varargout] = filterCatalog(tMain,ccMain,thresh,varargin)
arginmin = 3;
if nargin < arginmin
    fprintf(2,"not enough input arguments\n");
    return;
end

%%
[t,cc,removeIndices] = removeRepeatedMatches(tMain,ccMain,thresh);

inArgs = varargin;
varargout = inArgs;
if isempty(removeIndices)
    varargout = inArgs;
    return;
end

if nargin > arginmin
    for j = 1:nargin-arginmin
        extravar = inArgs{j};
        for i = 1:length(removeIndices)
            rI = removeIndices{i};
            extravar(rI,:) = [];
        end
        varargout{j} = extravar;
    end
end