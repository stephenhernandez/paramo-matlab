function [ccnorm,t] = updateCustomCatalog(S,templateFileName,temporaryCatalogFile,longTermCatalogFile,...
    threshold,maxTemplates,recordLength,maxN,diffFlag,mpd,...
    linearccnorm,plotFlag,verboseFlag)
% function [ccnorm,t] = updateCustomCatalog(S,templateFileName,temporaryCatalogFile,longTermCatalogFile,...
%     threshold,maxTemplates,startIndex,recordLength,maxN,diffFlag,mpd,subspaceFlag,...
%     linearccnorm,plotFlag,verboseFlag)


%% gen new
[~,tabs,NCC,z2p,Neff,p2rms,kurt,ccnorm,t] = subspaceDetector(S,...
    threshold,templateFileName,maxTemplates,recordLength,maxN,mpd,...
    diffFlag,linearccnorm,plotFlag,verboseFlag);

%%
fprintf("saving TEMPORARY file: %s\n",temporaryCatalogFile);
save(temporaryCatalogFile,'tabs','p2rms','z2p','NCC','Neff','kurt');
clear tabs p2rms z2p NCC Neff kurt

%% sync datasets
fprintf("loading: %s\n",longTermCatalogFile);
if ~exist(longTermCatalogFile,'file')
    load(temporaryCatalogFile,'tabs','p2rms','z2p','NCC','Neff','kurt');
    save(longTermCatalogFile,'tabs','p2rms','z2p','NCC','kurt','Neff');
    load(longTermCatalogFile,'tabs','p2rms','z2p','NCC','kurt','Neff');
else
    load(longTermCatalogFile,'tabs','p2rms','z2p','NCC','kurt','Neff');
end

%%
T = load(longTermCatalogFile,'tabs');
lastT = min([min(pull(S,'ref')) + seconds(60) max(T.tabs)]);
clear T;

%% this are new (temporary) data to append to long-term catalog
tNew = load(temporaryCatalogFile,'tabs');
tNew = tNew.tabs;
p2rmsNew = load(temporaryCatalogFile,'p2rms');
p2rmsNew = p2rmsNew.p2rms;
z2pNew = load(temporaryCatalogFile,'z2p');
z2pNew = z2pNew.z2p;
NCCNew = load(temporaryCatalogFile,'NCC');
NCCNew = NCCNew.NCC;
kurtNew = load(temporaryCatalogFile,'kurt');
kurtNew = kurtNew.kurt;
neffNew = load(temporaryCatalogFile,'Neff');
neffNew = neffNew.Neff;

%% filter here
newGood = tNew > lastT;
tKeep = tabs <= lastT;

%% append here
tabs = [tabs(tKeep); tNew(newGood)];
p2rms = [p2rms(tKeep,:); p2rmsNew(newGood,:)];
z2p = [z2p(tKeep,:); z2pNew(newGood,:)];
NCC = [NCC(tKeep); NCCNew(newGood)];
kurt = [kurt(tKeep,:); kurtNew(newGood,:)];
Neff = [Neff(tKeep); neffNew(newGood)];

%%
clear tNew p2rmsNew z2pNew NCCNew kurtNew neffNew lastT;

%%
fprintf("saving LONG-TERM file: %s\n",longTermCatalogFile);
save(longTermCatalogFile,'tabs','p2rms','z2p','NCC','kurt','Neff');