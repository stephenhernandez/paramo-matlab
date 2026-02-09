function d = bitsign(x,n)
%
% d = bitsign(x,n)
% returns signed double value from unsigned n-bit number x.
% more efficient version of bitsplit(x,n,n)
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019

%%
d = double(bitand(x,oldbitcmp32(0,n))) - double(bitget(x,n)).*(2^n);