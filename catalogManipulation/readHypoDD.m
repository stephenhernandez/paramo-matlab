function [t,eqlat,eqlon,eqdepth,eqmag,id] = readHypoDD(regionName,fullFileName)
% Inputs:
% regionName = 'ecuador';

%% parse inputs, deal variables
if nargin < 1
    regionName = 'sierra_negra';
end

if nargin < 2
    fullFileName = false;
end

baseDir = fullfile("~","masa","relocation");
if ~fullFileName
    file_loc = fullfile(baseDir,regionName,"hypodd","hypoDD.reloc");
else
    file_loc = regionName;
end

%% read sc3 catalog and filter
unix(['grep -v "\*" ',char(file_loc),' > tmp.txt']);
unix(['mv tmp.txt ',char(file_loc)]);

if ~fullFileName
    mdat = load(file_loc);
    id = mdat(:,1);
    eqmag = mdat(:,17);
    eqlon = mdat(:,3);
    eqlat = mdat(:,2);
    eqdepth = mdat(:,4);
    t = datetime(mdat(:,11),mdat(:,12),mdat(:,13),mdat(:,14),mdat(:,15),mdat(:,16));
else
    T = readtable(file_loc);
    t = datetime(table2array(T(:,2:7)));
    id = T.Var1;
    eqlat = T.Var8;
    eqlon = T.Var9;
    eqdepth = T.Var10;
    eqmag = T.Var11;
end
[t,sI] = sort(t);
eqlon = eqlon(sI);
eqlat = eqlat(sI);
eqdepth = eqdepth(sI);
eqmag = eqmag(sI);
id = id(sI);