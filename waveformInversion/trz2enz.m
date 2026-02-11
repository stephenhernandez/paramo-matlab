% function trz2enz(S,baz)
% baz: azimuth from receiver to source, in degrees
% S should have 2 components, first component is transverse, second component is radial
% components beyond the second are ignored

function S = trz2enz(S,baz)
if nargin < 2; baz = 0; end

%%
lS = length(S);
if lS < 2
    disp('struct must have at least two components');
    return;
end

%%
if baz >= 360
    baz = baz - 360;
end

if baz <= 180
    rotang = baz + 180;
else
    rotang = baz - 180;
end

% rotang is essentially the az (assuming locs not too close to any of the poles)
rotang = -rotang;

disp(rotang);
S = rotateSacData(S,rotang,1,2,false);

%%
kcmpnm = S(cmp1).kcmpnm;
kcmpnm = char(kcmpnm);
kcmpnm(3) = "E";
S(cmp1).kcmpnm = string(kcmpnm);
kcmpnm = S(cmp2).kcmpnm;
kcmpnm = char(kcmpnm);
kcmpnm(3) = "N";
S(cmp2).kcmpnm = kcmpnm;
