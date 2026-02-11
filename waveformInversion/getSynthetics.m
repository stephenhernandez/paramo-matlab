function S = getSynthetics(strike,dip,rake,evla,evlo,evdp,stla,stlo,scalar,prefix)
% fixed strike, dip, rake, and fixed evla,evlo
% allows multiple sensors (multiple distances,azimuths)

if nargin < 9; scalar = 1e6; end
if nargin < 10; prefix = '~/gf/ecuador/agudelo1Hz_triangleSTF_L2_Disp/'; end

%%
refEllipse = referenceEllipsoid('wgs84');
[dist,azs] = distance(evla,evlo,stla,stlo,refEllipse);
dist = dist*1e-3;

%%
dI = dist > 600;
if sum(dI)
    disp('some sensors are too far');
    dist = dist(~dI);
end

%%
ld = length(dist);
S = populateSacStruct(3);
S = repmat(S,1,ld);
for i = 1:ld
    S_ = forwardModel(strike,dip,rake,dist(i),evdp,azs(i),scalar,prefix);
    S(:,i) = S_;
end