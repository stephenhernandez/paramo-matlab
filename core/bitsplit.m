function d = bitsplit(x,b,n)
%
% d = bitsplit(x,b,n)
% splits the b-bit number x into signed n-bit array
% x: unsigned integer class
% n: ranges from 1 to b
% b: a multiple of n
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019

%%
rx = size(x,1);
sign = repmat((b:-n:n)',1,rx);
kshifts = repmat((n:n:b)' - b,1,rx);
x = repmat(x',b/n,1);
d = double(bitand(bitshift(x,kshifts),oldbitcmp32(0,n))) - double(bitget(x,sign))*(2^n);