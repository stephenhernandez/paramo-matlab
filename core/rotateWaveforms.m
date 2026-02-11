function S = rotateWaveforms(S,theta,eci,nci,RENAMEFLAG)
%
% S = rotateWaveforms(S,theta,eci,nci,RENAMEFLAG)
% rotates waveforms clockwise from north by theta degrees.
% need at least two components
%
% theta: angle in degrees [default 0 degrees]
% eci: index to easting channel [default 1]
% nci: index to northing channel [default 2]
% RENAMEFLAG: do i rename the components to R,T? [default false]
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
if nargin < 2
    theta = 0;
end

if nargin < 3
    eci = 1;
end

if nargin < 4
    nci = 2;
end

if nargin < 5
    RENAMEFLAG = false;
end

%%
Se = S(eci);
Sn = S(nci);
estart = Se.ref + Se.b; 
nstart = Sn.ref + Sn.b;
if estart ~= nstart
    fprintf("east and north component are not synced\n")
    return;
end

%%
e = Se.d;
pts_e = S(eci).npts;
n = Sn.d;
pts_n = S(nci).npts;

if pts_e > pts_n
    e = e(1:pts_n);     % cut e to be size of n
elseif pts_n > pts_e
    n = n(1:pts_e);     % cut n to be size of e
end
[erot,nrot] = rotate2d(e,n,theta);
S(eci) = dealHeader(Se,erot);
S(nci) = dealHeader(Sn,nrot);

%%
if ~RENAMEFLAG
    return;
end

kcmpnm = S(eci).kcmpnm;
kcmpnm = char(kcmpnm);
kcmpnm(3) = 'T';
S(eci).kcmpnm = string(kcmpnm);

kcmpnm = S(nci).kcmpnm;
kcmpnm = char(kcmpnm);
kcmpnm(3) = 'R';
S(nci).kcmpnm = string(kcmpnm);