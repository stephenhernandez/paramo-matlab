function d = taper(d,tw)
%
% taper returns waveforms with edges tapered with quarter wavelength of a cosine
%
% d = taper(d,tw)
% tw: taper width, expressed as decimal percentage (50% = 0.50, 25% on left side, 25% on right)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    tw = false;
end

if ~tw
    return;
end

%% for versions of matlab after 2016b, implicit expansion is done
r = size(d,1);
if tw > 1
    tw = tw/r;
end
scalar = tukeywin(r,tw);
d = d.*scalar;