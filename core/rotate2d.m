function [erot,nrot,RM] = rotate2d(echan,nchan,theta)
%
% rotate2d rotate waveformsto theta. need at least two components
%
% [erot,nrot] = rotate2d(echan,nchan,theta)
% theta: angle in degrees clockwise from north
% tested: 28 JAN 2026
% if rm = [costheta -sintheta; sintheta costheta]
% then [echan,nchan]*rm ==> erot = costheta*echan+sintheta*nchan 
% and                       nrot = -sintheta*echan+costheta*nchan
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
RM = rotation_matrix(theta);
R = [echan,nchan]*RM;
erot = R(:,1);
nrot = R(:,2);