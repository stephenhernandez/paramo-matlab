function E = readSCBulletin(events,diasFlag)
if nargin < 2
    diasFlag = false;
end

% 2017 11 17 13 40 41.023   0.107 -2.41423   0.989 -80.00224   1.136 61.94981   1.895   1.181  43.619 167 igepn2017wnrz locsat iasp91 manual confirmed igepn 6.237854893  0.19 42 mlv trimmedmean unk 2017 11 22 21 29 02
% 167 Phase arrivals:
% GYE3 EC HNE  0.278 21.8 S 2017 11 17 13 40 58.833 -2.646 M    1.00 0
% GYE3 EC HNZ  0.278 21.8 P 2017 11 17 13 40 50.460 -2.192 M    1.00 -1
% AGYE EC HNZ  0.361  8.0 P 2017 11 17 13 40 51.777 -1.628 M    1.00 -1
% APLA EC HNE  0.448 240.5 S 2017 11 17 13 41 05.359  1.147 M    1.00 0
% APLA EC HNZ  0.448 240.5 P 2017 11 17 13 40 53.484 -0.705 M    1.00 0
% MILO EC HHE  0.501 63.0 S 2017 11 17 13 41 03.748 -1.315 M    1.00 0
% MILO EC HHZ  0.501 63.0 P 2017 11 17 13 40 53.122 -1.546 M    1.00 -1
% AMIL EC HNE  0.527 63.9 S 2017 11 17 13 41 03.884 -1.595 M    1.00 0
% AMIL EC HNZ  0.527 63.9 P 2017 11 17 13 40 53.190 -1.712 M    1.00 -1
%
% 167 Station magnitudes:
% MILO EC BHZ  0.501 63.0 MLv  6.45  0.21 5723.09
% MILO EC BHZ  0.501 63.0 Mjma  6.08  0.08 2838.66
% ISPG EC BHZ  0.576 196.9 MLv  5.40 -0.84 377.588
% ISPG EC BHZ  0.576 196.9 Mjma  5.23 -0.78 311.432
% COHC EC BHZ  0.746 94.0 MLv  6.44  0.20 3339.4
% COHC EC BHZ  0.746 94.0 Mjma  7.04  1.03 12842.6
% SALI EC BHZ  1.014 282.9 MLv  5.84 -0.40 596.319
% SALI EC BHZ  1.014 282.9 Mjma  5.83 -0.18 468.318

% function E = readSCBulletin(event_file)

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Saturday, Jul 27, 2019

%%
if ~diasFlag
    lEvents = length(events);
    E = populateSeisCompStructure(lEvents);
    for i = 1:lEvents
        event_ = events(i);
        event_ = char(event_);
        if strcmp(event_(end-3:end),'.txt')
            event_ = event_(1:end-4);
        end
        yyyyStr = event_(end-7:end-4);
        event_ = strcat('~/phaseInformationSC5/',yyyyStr,'/',event_,'.txt');

        [t,lat,lon,depth,mag,timerr,laterr,lonerr,deptherr,magerr,rms,azgap,usedPhases,nmag,...
            id,methodID,earthModel,evMode,evStatus,agencyID,magtype,magmethod,evDescription,...
            creationTime,authorID,evType] = readSourceParameters(event_);
        [Pphases,nPphases] = readPhaseInformation(event_,true,t);
        [Sphases,nSphases] = readPhaseInformation(event_,false,t);
        [MLv,nMLv,mB,nmB,ML,nML,Mjma,nMjma,MsBB,nMsBB,Mwp,nMwp,mb,nmb,magerr2] = readMagnitudeInformation(string(event_));


        E(i).t = t;
        E(i).lat = lat;
        E(i).lon = lon;
        E(i).depth = depth;
        E(i).mag = mag;
        E(i).magerr2 = magerr2;
        E(i).timerr = timerr;
        E(i).laterr = laterr;
        E(i).lonerr = lonerr;
        E(i).deptherr = deptherr;
        E(i).magerr = magerr;
        E(i).rms = rms;
        E(i).azgap = azgap;
        E(i).usedPhases = usedPhases;
        E(i).nmag = nmag;

        E(i).id = id;
        E(i).methodID = methodID;
        E(i).earthModel = earthModel;
        E(i).evMode = evMode;
        E(i).evStatus = evStatus;
        E(i).agencyID = agencyID;
        E(i).magtype = magtype;
        E(i).magmethod = magmethod;
        E(i).evDescription = evDescription;
        E(i).authorID = authorID;
        E(i).evType = evType;

        E(i).creationTime = creationTime;
        E(i).nPphases = nPphases;
        E(i).Pphases = Pphases;
        E(i).nSphases = nSphases;
        E(i).Sphases = Sphases;
        E(i).nMLv = nMLv;
        E(i).MLv = MLv;
        E(i).nMjma = nMjma;
        E(i).Mjma = Mjma;
        E(i).nML = nML;
        E(i).ML = ML;
        E(i).nMsBB = nMsBB;
        E(i).MsBB = MsBB;
        E(i).nMwp = nMwp;
        E(i).Mwp = Mwp;
        E(i).nmB = nmB;
        E(i).mB = mB;
        E(i).nmb = nmb;
        E(i).mb = mb;
        %E = struct2table(E,'AsArray',true);
    end
else
    % Event:
    %     Public ID              igepn2014what
    %     Description
    %       region name: Colombia-Ecuador Border Region
    % Origin:
    %     Date                   2014-11-13
    %     Time                   22:07:32.653  +/- 4.925 s
    %     Latitude               0.7817 deg  +/-   38 km
    %     Longitude             -77.9325 deg  +/-    5 km
    %     Depth                    3.00 km   (fixed)
    %     Agency                 IGEPN
    %     Mode                   manual
    %     Status                 confirmed
    %     Residual RMS              0.1 s
    %     Azimuthal gap             270 deg
    %
    % 2 Network magnitudes:
    %     MLv       0.92 +/- 0.31   4
    %     M         0.92            4 preferred
    %
    % 5 Phase arrivals:
    %     sta  net   dist azi  phase   time         res     wt  sta
    %     CHL1  EC    0.0 307  P       22:07:33.561   0.0 M  1.0  CHL1
    %     CHL2  EC    0.0  35  P       22:07:33.586  -0.1 M  1.0  CHL2
    %     ECEN  EC    0.1 314  P       22:07:34.248  -0.0 M  1.0  ECEN
    %     CERN  OP    0.1 328  P       22:07:34.378   0.0 M  1.0  CERN
    %     IPAN  OP    0.1  37  P       22:07:34.934   0.1 M  1.0  IPAN
    %
    % 4 Station magnitudes:
    %     sta  net   dist azi  type   value   res   amp per
    %     CHL1  EC    0.0 307  MLv     1.44  0.52   1.2323
    %     CHL2  EC    0.0  35  MLv     0.70 -0.22   0.212206
    %     CERN  OP    0.1 328  MLv     0.94  0.02   0.297251
    %     IPAN  OP    0.1  37  MLv     0.75 -0.17   0.162896
    % function E = readSCBulletin(event_file)

    events = char(events);
    %     idBull = regexp(events,'.bulletin','once');
    %     if isempty(idBull)
    %         events = strcat('~/research/now/sierra_negra/NLL_TEPP_20190904/BULLETINS/',events,'.bulletin');
    %     else
    %         events = strcat('~/research/now/sierra_negra/NLL_TEPP_20190904/BULLETINS/',events);
    %     end
    [origyyyy,origmm,origdd,orighh,origmmm,origsec,origlat,origlaterr,origlon,origlonerr,origdepth,origrms,origgap,id] = ...
        readSourceParametersOld(events);
    [origmag1,nPphases,Pphases] = readPhaseInformationOld(events,origyyyy,origmm,origdd,orighh,origmmm,origsec);
    [origmag2,nSphases,Sphases] = readSphaseInformationOld(events,origyyyy,origmm,origdd,orighh,origmmm,origsec);
    [origmag3,nMLv,MLv,nML,ML,nMjma,Mjma,nMsBB,MsBB,nMwp,Mwp] = readMagnitudeInformationOld(events);

    %%
    lP = length(Pphases);
    for i = 1:lP
        stnm_ = char(Pphases(i).stnm);
        if strcmp(stnm_(1:2),'SN')
            Pphases(i).ntwk = "9D";
            Pphases(i).chan = "HHZ";
        elseif strcmp(stnm_,'CEAZ') || strcmp(stnm_,'FER1')
            Pphases(i).chan = "BHZ";
        else
            Pphases(i).chan = "HNZ"; %for cases like VCH1, PVIL, etc
        end
    end

    %%
    origmags = [origmag1 origmag2 origmag3];
    origmagI = isfinite(origmags);
    E = populateSeisCompStructure();
    E.mag = max(origmags(origmagI));
    E.nPphases = nPphases;
    E.Pphases = Pphases;
    E.nSphases = nSphases;
    E.Sphases = Sphases;
    E.nMLv = nMLv;
    E.MLv = MLv;
    E.nMjma = nMjma;
    E.Mjma = Mjma;
    E.nML = nML;
    E.ML = ML;
    E.nMsBB = nMsBB;
    E.MsBB = MsBB;
    E.nMwp = nMwp;
    E.Mwp = Mwp;
    E.lat =origlat;
    E.laterr =origlaterr;
    E.lon =origlon;
    E.lonerr = origlonerr;
    E.depth = origdepth;
    E.rms = origrms;
    E.azgap = origgap;
    E.t = datetime([origyyyy origmm origdd orighh origmmm origsec]);
    E.id = string(id);
end