function [center_lat,center_lon,minLat,maxLat,minLon,maxLon,ewindex,nsindex]...
    = get_region_dimensions(regionName)
minLat=[];
maxLat=[];
minLon=[];
maxLon=[];
center_lon = [];
center_lat = [];
ewindex = [];
nsindex = [];

if strcmp(regionName,'cotopaxi')
    minLat=-0.82;
    maxLat=-0.55;
    minLon=-78.57;
    maxLon=-78.30;
    center_lon = -78.4361;%-78.2165
    center_lat = -0.6836;
    ewindex = 493;
    nsindex = 495;
elseif strcmp(regionName,'pichincha') || strcmp(regionName,'pichincha_dem') || strcmp(regionName,'pichincha_30m_dem')
%     minLat = -0.26; %-0.32;
%     maxLat = -0.09; %-0.01;
%     minLon = -78.66; %-78.73;
%     maxLon = -78.54; %-78.45;
    minLat = -0.24;
    maxLat = -0.1;
    minLon = -78.70;
    maxLon = -78.52;
    center_lon = -78.6124;
    center_lat = -0.1705;
    ewindex = 591;
    nsindex = 436;
elseif strcmp(regionName,'wolf')
    minLat = -0.26; %-0.32;
    maxLat = -0.09; %-0.01;
    minLon = -78.66; %-78.73;
    maxLon = -78.54; %-78.45;
    center_lon = -91.343477;
    center_lat = 0.01936;
    ewindex = 591;
    nsindex = 436;    
elseif strcmp(regionName,'reventador')
    minLat=-0.18;
    maxLat=0.01;
    minLon=-77.75;
    maxLon=-77.55;
    center_lon = -77.6581;
    center_lat = -0.0810;
    ewindex = 336;
    nsindex = 341;
elseif strcmp(regionName,'esmeraldas')
    %-79.888504                -79.358513                  0.604521                  1.129349
    minLat = 0.58;
    maxLat = 1.2;
    minLon = -80.02;
    maxLon = -79.3;
    center_lon = -77.6581;
    center_lat = -0.0810;
    ewindex = 336;
    nsindex = 341;
elseif strcmp(regionName,'atacames')
    minLat = 0.45;
    maxLat = 1.2;
    minLon = -80.8;
    maxLon = -79.7;
    center_lon = -77.6581;
    center_lat = -0.0810;
    ewindex = 336;
    nsindex = 341;
elseif strcmp(regionName,'puntaGaleras')
    minLat = 0.;
    maxLat = 1.2;
    minLon = -80.8;
    maxLon = -79.5;
    center_lon = -77.6581;
    center_lat = -0.0810;
    ewindex = 336;
    nsindex = 341;
elseif strcmp(regionName,'tungurahua')
    minLat=-1.54;
    maxLat=-1.35;
    minLon=-78.54;
    maxLon=-78.35;
    center_lon = -78.445;
    center_lat = -1.4702;
    ewindex = 434;
    nsindex = 340;
elseif strcmp(regionName,'cuicocha')
    minLat=0.20; %0.24;
    maxLat=0.45; %0.46;
    minLon=-78.4; %-78.45;
    maxLon=-78.2; %-78.21;
    center_lon = -78.3479;
    center_lat = 0.3610;
    ewindex = 364;
    nsindex = 378;
elseif strcmp(regionName,'atuntaqui')
    minLat=0.20; %0.24;
    maxLat=0.45; %0.46;
    minLon=-78.4; %-78.45;
    maxLon=-78.2; %-78.21;
    center_lon = -78.3479;
    center_lat = 0.3610;
    ewindex = 364;
    nsindex = 378;    
elseif strcmp(regionName,'chiles')
    minLat=0.74; %0.6;
    maxLat=0.84; %0.9;
    minLon=-78.0;
    maxLon=-77.9; %-77.75;
    %minLat=0.74;
    %maxLat=0.91;
    %minLon=-78.02;
    %maxLon=-77.88;
    center_lon = -77.938;
    center_lat = 0.8170;
    ewindex = 297;
    nsindex = 336;
elseif strcmp(regionName,'chilesExpanded')
    minLat=0.65; %0.6; 
    maxLat=0.9;
    minLon=-78.0;
    maxLon=-77.75;    
    center_lon = -77.938;
    center_lat = 0.8170;
    ewindex = 297;
    nsindex = 336;
elseif strcmp(regionName,'sangay')
    minLat=-2.11;
    maxLat=-1.9;
    minLon=-78.45;
    maxLon=-78.23;
    center_lon = -78.3421;
    center_lat = -2.0051;
    ewindex = 388;
    nsindex = 400;
elseif strcmp(regionName,'pedernales')
    minLon = -82;
    maxLon = -79;
    minLat = -2.;
    maxLat = 2;
    %     minLat=-6;
    %     maxLat=4;
    %     minLon=-86;
    %     maxLon=-75;
    center_lon = -78.0;
    center_lat = -2.0;
    ewindex = 0;
    nsindex = 0;
elseif strcmp(regionName,'cayambe')
    minLat=-0.07;
    maxLat=0.15;
    minLon=-78.15;
    maxLon=-77.85;
    center_lon = -77.9872;
    center_lat = 0.025;
    ewindex = 534;
    nsindex = 604;
elseif strcmp(regionName,'quito')
    minLat=-0.3;
    maxLat=-0.1; %0.0;
    minLon=-78.55; %-78.7;
    maxLon=-78.4; %-78.2;
    center_lon = -77.9872;
    center_lat = 0.025;
    ewindex = 534;
    nsindex = 604;
elseif strcmp(regionName,'quito_valley')
    minLat=-0.25;
    maxLat=0.0;
    minLon=-78.52;
    maxLon=-78.32;
    center_lon = -77.9872;
    center_lat = 0.025;
    ewindex = 534;
    nsindex = 604;
elseif strcmp(regionName,'galapagos')
    minLat=-2;
    maxLat=1.0;
    minLon=-92;
    maxLon=-89;
    center_lon = -91.39;
    center_lat = -0.928;
    ewindex = 2045;
    nsindex = 739;
elseif strcmp(regionName,'southernIsabela')
    minLat=-1.1;
    maxLat=-0.6;
    minLon=-91.35;
    maxLon=-90.9;
    center_lon = -91.131998;
    center_lat = -0.814308;
    ewindex = 2045;
    nsindex = 739;
elseif strcmp(regionName,'sierraNegra')
    minLat=-1.2;
    maxLat=-0.;
    minLon=-91.85;
    maxLon=-90.75;
    center_lon = -91.131998;
    center_lat = -0.814308;
    ewindex = 2045;
    nsindex = 739;
elseif strcmp(regionName,'bahiaElizabeth')
    minLat = -0.725; %-1;
    maxLat = -0.5; %0 ; %-0.6;
    minLon = -91.35; %-91.3;
    maxLon = -91.05; %-90.95;
    center_lon = -91.131998;
    center_lat = -0.814308;
    ewindex = 2045;
    nsindex = 739;
elseif strcmp(regionName,'sierraNegraCaldera') || strcmp(regionName,'sierra_negra')
    minLat = -1.1; %-0.865; %-0.95;
    maxLat = -0.5; %-0.765; %-0.65;
    minLon = -91.5; %-91.2; %-91.25;
    maxLon = -90.9; %-91.07; %-91;
    center_lon = -91.131998;
    center_lat = -0.814308;
    ewindex = 2045;
    nsindex = 739;    
elseif strcmp(regionName,'cerro_azul')
    minLat=-1.25;
    maxLat=-0.25;
    minLon=-92;
    maxLon=-90.5;
    center_lon = -91.39;
    center_lat = -0.928;
    ewindex = 2045;
    nsindex = 739;
elseif strcmp(regionName,'fernandina')
    minLat=-0.65; 
    maxLat=-0.05;
    minLon=-91.8;
    maxLon=-91.2;
    center_lon = -91.54;
    center_lat = -0.372;
    ewindex = 2045;
    nsindex = 739;
elseif strcmp(regionName,'balao') || strcmp(regionName,'balao_dem')
    minLon = -80; %80.15;
    maxLon = -79.5; %45;
    minLat = -3.2; %-3.25;
    maxLat = -2.5; %35; %balao
    center_lon = -79.8;
    center_lat = -2.9;
    ewindex = 2045;
    nsindex = 739;
elseif strcmp(regionName,'balaoCutDEM')
    minLon = -80.3;
    maxLon = -79.3;
    minLat = -3.3;
    maxLat = -2.3; %balao
    center_lon = -79.8;
    center_lat = -2.9;
    ewindex = 2045;
    nsindex = 739;
else
    disp('sorry, no region with that name found.')
end