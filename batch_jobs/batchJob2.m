% clear; close all; clc;
% [eqType,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,nPhases,nMLv,...
% timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
% locMethod,earthModel,creationTime,agencyID,evMode,scbullmag] = readCat1(datetime(2010,01,01),datetime(2050,01,01),0);
% boundaryBox = getRegionSpatialDimensions('ecuador');
% err = sqrt(eqdeptherr.^2 + eqlonerr.^2 + eqlaterr.^2);
% err2 = sqrt(timerr.^2 + stderr.^2);
% eI = nPhases >= 4 & azgap <= 225 & err <= 10 & err2 < 2;
%
%
% close all;
% plotEcuadorSeismicity(t(eI),[eqlat(eI) eqlon(eI) eqdepth(eI) eqmag(eI)],boundaryBox,true,true);

%%
clear; close all; clc;
[eqType,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,nPhases,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,creationTime,agencyID,evMode,scbullmag] = readCat1();
tI = t <= datetime(2020,10,02);
tI = find(tI);
for i = length(tI):-1:1
    loadEventWaveforms(ids(tI(i)));
end