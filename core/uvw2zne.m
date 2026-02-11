function S = uvw2zne(S,make)
if nargin < 2
    make = "reftek";
end

if strcmp(make,"reftek")    %colt
    RM = [-2 1 1;...
        0 sqrt(3) -sqrt(3);...
        sqrt(2) sqrt(2) sqrt(2)];
else                        %trillium?
    RM = [2 -1 -1;...
        0 sqrt(3) -sqrt(3);...
        sqrt(2) sqrt(2) sqrt(2)];
end
RM = RM/sqrt(6);
d = pull(S);
d = d';
xyz = RM*d;
xyz = xyz';

S(1) = dealHeader(S(1),xyz(:,2));   % E
S(2) = dealHeader(S(2),xyz(:,3));   % N
S(3) = dealHeader(S(3),xyz(:,1));   % Z
clear d xyz;

newChans = ["E";"N";"Z"];
for i = 1:3
    kcmpnm_ = S(i).kcmpnm;
    kcmpnm_ = char(kcmpnm_);
    kcmpnm_ = kcmpnm_(1:2);
    kcmpnm_ = strcat(string(kcmpnm_),newChans(i));
    S(i).kcmpnm = kcmpnm_;
end
S = flipud(S);