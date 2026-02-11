function S = populateRSAM(lfcs,hfcs,meanFlags,rmsFlags)
% if nargin < 1
%     sizeS = 1;
% end

%%
S.t = NaN;

%%
techniques = rsamTechniqueStr(meanFlags,rmsFlags);
techniques = unique(techniques);
ltechniques = length(techniques);

llfcs = length(lfcs);
for i = 1:llfcs
    lfc = lfcs(i);
    hfc = hfcs(i);

    for j = 1:ltechniques
        technique = techniques(j);
        fieldName = rsamFieldNames(lfc,hfc,technique);
        S.(fieldName) = NaN;
    end
end

%%
% if numel(sizeS) == 1
%     S = repmat(S,sizeS,1); %make column vector
% else
%     S = repmat(S,sizeS);
% end
