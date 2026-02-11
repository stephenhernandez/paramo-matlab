function output = oldbitcmp(x,n)
%
% output = oldbitcmp(x,n)
%
% oldbitcmp mimics the behavior of bitcmp versions of yore
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019

%%
if nargin < 2
    output = bitcmp(x);
    return;
end

%%
maxN = 2^n-1;   %max number one can represent in n-bits
fmtCode = 0;
switch fmtCode
    case 0 % uint32
        fmt  = 'uint32';
        out1 = bitcmp(uint32(x),fmt);
        out2 = bitcmp(uint32(maxN),fmt);
        output = out1 - out2;
    case 1 % uint8
        fmt  = 'uint8';
        out1 = bitcmp(uint8(x),fmt);
        out2 = bitcmp(uint8(maxN),fmt);
        output = out1 - out2;
    case 2 % uint16
        fmt  = 'uint16';
        out1 = bitcmp(uint16(x),fmt);
        out2 = bitcmp(uint16(maxN),fmt);
        output = out1 - out2;
end