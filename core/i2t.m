function t = i2t(i,ref,delta)
%
% i2t index-to-time
% t = i2t(i,ref,delta)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
t = ref + seconds((i - 1)*delta);