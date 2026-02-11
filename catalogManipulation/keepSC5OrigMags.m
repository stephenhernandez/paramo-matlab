function [sc5mag,sc5id,scbullmag,scbullid] = keepSC5OrigMags()
sc5orig = importdata('~/phaseInformationSC5/sc5_mags.txt');     % this is a preferred mag
scbull = importdata('~/phaseInformationSC5/scbull_mags.txt');  % this is an MLv

%% just the unique ids/magsin both sets
sc5id = string(sc5orig.textdata);
sc5mag = sc5orig.data;
[sc5id,ia] = unique(sc5id);
sc5mag = sc5mag(ia);

%%
badI = sc5mag == -999;
sc5mag(badI) = [];
sc5id(badI) = [];

%%
scbullid = string(scbull.textdata);
scbullmag = scbull.data;
[scbullid,ia] = unique(scbullid);
scbullmag = scbullmag(ia);

%% sync the two catalogs
[lia,locb] = ismember(sc5id,scbullid);
sc5id = sc5id(lia);
sc5mag = sc5mag(lia);
scbullid = scbullid(locb(lia));
scbullmag = scbullmag(locb(lia));