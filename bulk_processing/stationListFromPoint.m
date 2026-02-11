function [kstnm,d,stla,stlo,stel,sI,SNCLs] = stationListFromPoint(stla_,stlo_)

% stationFile = '~/igdata/TODAS_ESTACIONES_TODAS.txt';
% fileID = fopen(stationFile);
% C1 = textscan(fileID,'%s %f %f %f');
% fclose(fileID);
% [eckstnm,ecstla,ecstlo,ecstel] = deal(C1{:});
% 
% % ec_meta_data = load('~/igdata/ecuador_meta_data_3');
% % eckstnm = ec_meta_data.kstnm;
% % ecstla = ec_meta_data.stla;
% % ecstlo = ec_meta_data.stlo;
% % ecstel = ec_meta_data.stel; %elevation in meters

%%
load('~/igdata/ecuadorSensorDataTable10','kstnm','knetwk','kcmpnm','khole','stla','stlo','stel','allSNCLs');

%%
refEllipse = referenceEllipsoid('wgs84');
d = distance(stla_,stlo_,stla,stlo,refEllipse)*1e-3;
[d,sI] = sort(d);
kstnm = string(kstnm(sI));
stla = stla(sI);
stlo = stlo(sI);
stel = stel(sI);
SNCLs = allSNCLs(sI);