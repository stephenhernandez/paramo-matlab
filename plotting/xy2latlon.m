function [lonRecovered,latRecovered] = xy2latlon(x,y,reflon,reflat,refEllipse)
if nargin < 5; refEllipse = referenceEllipsoid('wgs84'); end

refDist = x.^2 + y.^2;
refDist = sqrt(refDist);
refAz = atan2d(y,x);
refAz = 90-refAz;
[latRecovered,lonRecovered] = reckon(reflat,reflon,refDist*1e3,refAz,refEllipse);