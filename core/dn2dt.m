function tdt = dn2dt(tdn)
%
% dn2dt convert vector in datenum format to datetime object
%
% tdt = dn2dt(tdn)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
tdt = datetime(tdn,'ConvertFrom','datenum');