function y = oldbitcmp32(x,n)
%
% output = oldbitcmp32(x,n)
%
% oldbitcmp32 mimics the behavior of bitcmp versions of yore
%
% this is a stripped down version that specifically assumes unsigned
% 32-bit input
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019

%%
maxN = 2^n-1;   %max number one can represent in n-bits
y1 = bitcmp(uint32(x),'uint32');
y2 = bitcmp(uint32(maxN),'uint32');
y = y1 - y2;