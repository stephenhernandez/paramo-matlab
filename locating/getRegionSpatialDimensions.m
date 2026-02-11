function boundaryBox = getRegionSpatialDimensions(regionName)
minLat = [];
maxLat = [];
minLon = [];
maxLon = [];

if strcmp(regionName,'cotopaxi')
    minLat=0.65;
    maxLat=0.85;
    minLon=-78;
    maxLon=-77.77;

    % minLat=-1.6; %-0.78;
    % maxLat= 0.6; %-0.60;
    % minLon=-80;
    % maxLon=-76;
elseif strcmp(regionName,'lamerced')
    minLat=-0.65; %65456;
    maxLat=-0.55; %-0.5; %5;
    minLon=-78.55; %-78.5; %45; %3438;
    maxLon=-78.3; %34438;
elseif strcmp(regionName,'ecuador-continental') || strcmp(regionName,'ecuador')
    minLat = -6;
    maxLat = 2; %5.25;
    minLon = -83;
    maxLon = -74;
elseif strcmp(regionName,"seamount_manus")
    minLat = 0.16;
    maxLat = 0.4895;
    minLon = -80.2;
    maxLon = -79.80;
elseif strcmp(regionName,"napo")
    minLat = -1; %-0.9778;
    maxLat = -0.65; %-0.7778;
    minLon = -78.3; %-78.21661;
    maxLon= -77.8; %-78.01661;
elseif strcmp(regionName,'chingual')
    minLat = 0.3;
    maxLat = 0.75;
    minLon = -77.75;
    maxLon = -77.3;
elseif strcmp(regionName,'guagua') || strcmp(regionName,'pichincha')
    minLat = -0.24;
    maxLat = -0.1; %-0.09;
    minLon = -78.70;
    maxLon = -78.5;
elseif strcmp(regionName,'sumaco')
    centerLat = -0.541249;
    centerLon = -77.627416;
    dh = 0.25;
    minLat = centerLat - dh;
    maxLat = centerLat + dh;
    minLon = centerLon  - dh;
    maxLon = centerLon + dh;
elseif strcmp(regionName,'antisana')
    minLat = -0.65;
    maxLat = -0.35;
    minLon = -78.32;
    maxLon = -77.96;
elseif strcmp(regionName,'manta')
    minLat = -1.25;
    maxLat = -0.65;
    minLon = -81.05;
    maxLon = -80.35;
elseif strcmp(regionName,'plata')
    minLat = -1.75;
    maxLat = -1.1;
    minLon = -81.9;
    maxLon = -80.1;
elseif strcmp(regionName,'jama')
    minLat = -1.4;
    maxLat = 0.1;
    minLon = -81.3;
    maxLon = -80.1;
elseif strcmp(regionName,'outer-rise')
    minLat = -1.1;
    maxLat = 0.1;
    minLon = -82;
    maxLon = -80.8;
elseif strcmp(regionName,'wolf')
    minLat = -0.15;
    maxLat = 0.17;
    minLon = -91.6;
    maxLon = -91.15;
elseif strcmp(regionName,'darwin')
    minLat = -0.34;
    maxLat = -0.07;
    minLon = -91.45;
    maxLon = -91.10;
elseif strcmp(regionName,'alcedo')
    minLat = -0.6;
    maxLat = -0.3;
    minLon = -91.28;
    maxLon = -90.92;
elseif strcmp(regionName,'reventador')
    minLat= -0.16;
    maxLat= 0.0;
    minLon= -77.75;
    maxLon= -77.55;
elseif strcmp(regionName,'esmeraldas')
    minLat = 0.65; %0.45;
    maxLat = 1.45; %1.8;
    minLon = -80.2; %-80.4;
    maxLon = -79.1;
elseif strcmp(regionName,'atacames')
    minLat = 0.7;
    maxLat = 1.5;
    minLon = -80.5;
    maxLon = -79.7;
elseif strcmp(regionName,'punta-galera')
    minLat = 0.3;
    maxLat = 1.2;
    minLon = -80.9;
    maxLon = -79.6;
elseif strcmp(regionName,'tungurahua')
    minLat=-1.6;
    maxLat=-1.3;
    minLon=-78.6;
    maxLon=-78.3;
elseif strcmp(regionName,'pisayambo')
    minLat=-1.3;
    maxLat=-0.9;
    minLon=-78.5;
    maxLon=-78.0;
elseif strcmp(regionName,'chimborazo')
    minLat=-1.6;
    maxLat=-1.35;
    minLon=-78.98;
    maxLon=-78.66;
elseif strcmp(regionName,'quilotoa')
    maxLat=-0.78;
    minLat=-0.94;
    minLon=-79;
    maxLon=-78.8;
elseif strcmp(regionName,'cuicocha')
    minLat=0.24;
    maxLat=0.42;
    minLon=-78.42;
    maxLon=-78.22;
elseif strcmp(regionName,'atuntaqui')
    minLat = 0.20;
    maxLat = 0.45;
    minLon = -78.4;
    maxLon = -78.2;
elseif strcmp(regionName,'chiles')
    minLat=0.65;
    maxLat=0.85;
    minLon=-78;
    maxLon=-77.77;
elseif strcmp(regionName,'sangay')
    minLat= -2.08;
    maxLat= -1.92;
    minLon= -78.4;
    maxLon= -78.25;
elseif strcmp(regionName,'costa')
    minLon = -82;
    maxLon = -79;
    minLat = -2.;
    maxLat = 2;
elseif strcmp(regionName,'pululahua')
    minLon = -78.60;
    maxLon = -78.30;
    minLat = -0.11;
    maxLat = 0.13;
elseif strcmp(regionName,'cayambe')
    minLat = -0.06;
    maxLat = 0.18;
    minLon = -78.15;
    maxLon = -77.8;
elseif strcmp(regionName,'quito')
    minLat = -0.35;
    maxLat = 0.0;
    minLon = -78.55;
    maxLon = -78.25;
elseif strcmp(regionName,'puembo')
    minLat = -0.24; %-0.25;
    maxLat = -0.1; %-0.0;
    minLon = -78.46; %-78.55;
    maxLon = -78.30; %-78.30;
elseif strcmp(regionName,'quito_valley')
    minLat = -0.25;
    maxLat = 0.0;
    minLon = -78.52;
    maxLon = -78.32;
elseif strcmp(regionName,'galapagos')
    minLat = -1.8;
    maxLat = 1; %0.75;
    minLon = -92;
    maxLon = -88; %89;
elseif strcmp(regionName,'sierra_negra')
    minLat = -0.98;
    maxLat = -0.64;
    minLon = -91.3;
    maxLon = -90.92;
elseif strcmp(regionName,'bahia-elizabeth')
    minLat = -0.725;
    maxLat = -0.45;
    minLon = -91.45;
    maxLon = -91.05;
elseif strcmp(regionName,'cerro-azul') || strcmp(regionName,'cerro_azul')
    minLat = -1.10;
    maxLat = -0.75;
    minLon = -91.55;
    maxLon = -91.22;
elseif strcmp(regionName,'fernandina')
    minLat = -0.55;
    maxLat = -0.2;
    minLon = -91.7;
    maxLon = -91.3;
elseif strcmp(regionName,'balao')
    minLon = -80;
    maxLon = -79.5;
    minLat = -3.2;
    maxLat = -2.5;
elseif strcmp(regionName,'machala')
    minLon = -80.5;
    maxLon = -79.5;
    minLat = -3.5;
    maxLat = -3;
elseif strcmp(regionName,'puna')
    minLon = -80.7;
    maxLon = -79.3;
    minLat = -3.4;
    maxLat = -2.3;
elseif strcmp(regionName,'guayaquil')
    minLon = -80.35;
    maxLon = -79.4;
    minLat = -3;
    maxLat = -1.7;
elseif strcmp(regionName,'pedernales')
    minLat = -0.33;
    maxLat = 0.72;
    minLon = -81.25;
    maxLon = -79.55;
elseif strcmp(regionName,'machachi')
    minLat = -0.62;
    maxLat = -0.40;
    minLon = -78.56;
    maxLon = -78.38;
elseif strcmp(regionName,'pita')
    minLat = -0.62;
    maxLat = -0.40;
    minLon = -78.58;
    maxLon = -78.32;
elseif strcmp(regionName,'corazon')
    minLat = -0.62;
    maxLat = -0.4;
    minLon = -78.81;
    maxLon = -78.55;
else
    disp('sorry, no region with that name found.')
    boundaryBox = [minLon; maxLon; minLat; maxLat];
    return;
end

%%
boundaryBox = [minLon; maxLon; minLat; maxLat];