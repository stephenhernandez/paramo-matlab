% function S = forwardModel(strike,dip,rake,dist,depth,az,scalar,prefix)

% get synthetics given a strike,dip,rake
% output is in T,R,Z form; use trz2enz to switch

%%
function S = forwardModel(strike,dip,rake,dist,depth,az,scalar,db)
if nargin < 6; az = 0; end
if nargin < 7; scalar = 1e6; end
if nargin < 8; db = '~/gf/ecuador/agudelo1Hz_triangleSTF_L2_Disp/'; end

%% get all greens functions
[T,R,Z] = getGreensFunctionsAll(dist,depth,az,scalar,db);
T1 = pull(T);
R1 = pull(R);
Z1 = pull(Z);

%% get moment tensor components from strike dip and rake (uses garrets code)
mt = sdr2mt(strike,dip,rake);
mt(3) = []; %delete third (mpp) because is negative of sum of first two (zero trace)
mt = mt';

%%
T(1).d = T1*mt; 
T(1).('kcmpnm') = 'T'; 
T(1).('knetwk') = 'SY';

R(1).d = R1*mt; 
R(1).('kcmpnm') = 'R'; 
R(1).('knetwk') = 'SY';

Z(1).d = Z1*mt; 
Z(1).('kcmpnm') = 'Z'; 
Z(1).('knetwk') = 'SY';

%%
S = [T(1); R(1); Z(1)];
