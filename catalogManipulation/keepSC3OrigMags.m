function [sc5mag,sc5id] = keepSC3OrigMags()
sc5orig = importdata('~/phaseInformationSC5/sc5_mags.txt');     % this is a preferred mag
sc5bull = importdata('~/phaseInformationSC5/scbull_mags.txt');  % this is an MLv

%% just the unique ids/magsin both sets
sc5id = string(sc5orig.textdata);
sc5mag = sc5orig.data;
[sc5id,ia] = unique(sc5id);
sc5mag = sc5mag(ia);

scbullid = string(sc5bull.textdata);
scbullmag = sc5bull.data;
[scbullid,ia] = unique(scbullid);
scbullmag = scbullmag(ia);

%%
[lia,locb] = ismember(sc5id,scbullid);
sc5id = sc5id(lia);
sc5mag = sc5mag(lia);
scbullid = scbullid(locb(lia));
scbullmag = scbullmag(locb(lia));

%% check once
badI = sc5mag == -999;
sc5mag(badI) = [];
sc5id(badI) = [];

%%
lia = ismember(sc5id,scbullid);
sc5id = sc5id(lia);
sc5mag = sc5mag(lia);

lia2 = ismember(scbullid,sc5id);
scbullid = scbullid(lia2);
scbullmag = scbullmag(lia2);

% %% check again
% badI = sc3mag == -999;
% if sum(badI)
%     disp('rare case, both are bad');
%     disp(sc3id(badI))
%     sc3mag(badI) = [];
%     sc3id(badI) = [];
% end