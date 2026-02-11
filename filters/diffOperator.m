function iomega = diffOperator(n,Fs)
%
% diffOperator returns a n-point frequency-domain differentiator
%
% iomega = diffOperator(n,Fs)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 1
    n = 100;
end

if nargin < 2
    Fs = 1;
end

%%
f = Fs.*(0:n-1)'/n;
iomega = 1i*2*pi.*f;