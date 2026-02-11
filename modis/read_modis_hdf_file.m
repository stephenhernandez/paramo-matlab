function [lats,lons,t,qa,first_days,last_days,uncert] = read_modis_hdf_file(modis_hdf_file)

lats = [];
lons = [];
qa = [];
t = [];
first_days = t;
last_days = t;
uncert = [];

try
    burn_date = hdfread(modis_hdf_file, 'Burn Date');
catch
    return;
end

bI = burn_date>=1;
burned_pixels = sum(sum(bI));

if ~burned_pixels
    fprintf('no data for file: %s, number of burn pixels: %d\n',modis_hdf_file,burned_pixels);
    return;
end

%%
R = 6371007.181;
T = 1111950;
xmin = -20015109;
ymax = 10007555;
w = 463.31271653; %T/2400;

%%
fname_str = split(modis_hdf_file,".");
fname_str = fname_str{3};
H = str2double(fname_str(2:3));
V = str2double(fname_str(5:6));

burn_date_uncertainty = hdfread(modis_hdf_file,'Burn Date Uncertainty');
qa = hdfread(modis_hdf_file,'QA');
fd = hdfread(modis_hdf_file,'First Day');
ld = hdfread(modis_hdf_file,'Last Day');
yyyy = hdfread(modis_hdf_file,'year');
yyyy = yyyy{:};

[nRows,nCols] = size(burn_date);

i = (0:nRows-1)';
j = (0:nCols-1)';
[X,Y] = meshgrid(j,i);

yvec = ymax - (Y + 0.5)*w - V*T;
xvec = (X + 0.5)*w + H*T + xmin;

% wkt = "PROJCS[""MODIS Sinusoidal"",BASEGEOGCRS[""User"",DATUM[""World Geodetic Survey 1984"",SPHEROID[""Authalic_Spheroid"",6371007.181,0.0]],PRIMEM[""Greenwich"",0.0],UNIT[""Degree"",0.0174532925199433]],PROJECTION[""Sinusoidal""],PARAMETER[""False_Easting"",0.0],PARAMETER[""False_Northing"",0.0],PARAMETER[""Central_Meridian"",0.0],UNIT[""Meter"",1.0]]";
% p = projcrs(wkt);
% [latvec,lonvec] = projinv(p,xvec,yvec);

phi = yvec./R;
latvec = rad2deg(phi);
lonvec = rad2deg(xvec./(R*cos(phi)));

%%
lats = latvec(bI);
lons = lonvec(bI);
uncert = burn_date_uncertainty(bI);
tOrig = burn_date(bI);
t = datetime(yyyy,01,tOrig);
badI = isnat(t);
qa = qa(bI);
first_days = fd(bI);
last_days = ld(bI);
prevYear = first_days > tOrig;
nextYear = last_days < tOrig;
fd = datetime(yyyy,01,first_days);
ld = datetime(yyyy,01,last_days);
ld(nextYear) = datetime(yyyy+1,01,last_days(nextYear));
fd(prevYear) = datetime(yyyy-1,01,first_days(prevYear));
last_days = ld;
first_days = fd;

allLengths = [length(t); length(qa); length(uncert); length(first_days); length(last_days)];
if sum(badI)
    fprintf('difference_exists: %d, %s\n\n',sum(badI),modis_hdf_file);
end

