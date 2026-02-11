function [T,R,Z] = getGreensFunctions(dist,depth,az,scalar,prefix)
% fixed distance, fixed depth greens functions
% this is a non-volumetric case, where the trace is 0
if nargin < 3; az = 0; end
if nargin < 4; scalar = 1e6; end
if nargin < 5; prefix = '~/gf/ecuador/agudelo1Hz_triangleSTF_L2_Disp/'; end

%% nearest depth
if mod(depth,1) == 0
    depth = depth + 0.5;
end

%% find nearest distance (dist should be in km)
dist = round(dist)*10; %nearest distance
gstr = num2str(dist);
if dist < 100
    gstr = strcat('00',gstr);
elseif dist < 1000
    gstr= strcat('0',gstr);
end

%%
dpstr = num2str(depth);
if depth < 10
    dpstr = strcat('0',dpstr);
end
dpstr = strcat('H0',dpstr,'/');
base = strcat('GF.',gstr,'.SY.');

%% get raw synthetics
disp(strcat(prefix,dpstr,'RR/',base,'LHZ.SAC'));
S(1) = sac2struct(strcat(prefix,dpstr,'RR/',base,'LHZ.SAC'));
S(2) = sac2struct(strcat(prefix,dpstr,'TT/',base,'LHZ.SAC'));
S(3) = sac2struct(strcat(prefix,dpstr,'PP/',base,'LHZ.SAC'));
S(4) = sac2struct(strcat(prefix,dpstr,'RT/',base,'LHZ.SAC'));
S(5) = sac2struct(strcat(prefix,dpstr,'RR/',base,'LHL.SAC'));
S(6) = sac2struct(strcat(prefix,dpstr,'TT/',base,'LHL.SAC'));
S(7) = sac2struct(strcat(prefix,dpstr,'PP/',base,'LHL.SAC'));
S(8) = sac2struct(strcat(prefix,dpstr,'RT/',base,'LHL.SAC'));
S(9) = sac2struct(strcat(prefix,dpstr,'RP/',base,'LHT.SAC'));
S(10) = sac2struct(strcat(prefix,dpstr,'TP/',base,'LHT.SAC'));

%%
streams = pull(S);
uzrr = streams(:,1);
uztt = streams(:,2);
uzpp = streams(:,3);
uzrt = streams(:,4);
usrr = streams(:,5);
ustt = streams(:,6);
uspp = streams(:,7);
usrt = streams(:,8);
uerp = streams(:,9);
uetp = streams(:,10);

%%
Z = repmat(S(1),5,1);
R = repmat(S(5),5,1);
T = repmat(S(9),5,1);

%% rotate to azimuth
a = sin(az); b = cos(az);
c = sin(2*az); d = cos(2*az);
a2 = a*a;
b2 = b*b;
a2b2 = a2 - b2;

%% do the math
Z(1).d =  scalar*(uzrr - a2*uztt - b2*uzpp);
Z(2).d =  scalar*(d*uztt + a2b2*uzpp);
Z(3).d =  scalar*(b*uzrt);
Z(4).d = -scalar*(a*uzrt);
Z(5).d = -scalar*(c*(uztt-uzpp));

R(1).d =  scalar*(usrr - a2*ustt - b2*uspp);
R(2).d =  scalar*(d*ustt + a2b2*uspp);
R(3).d =  scalar*(b*usrt);
R(4).d = -scalar*(a*usrt);
R(5).d = -scalar*(c*(ustt-uspp));

T(1).d = -scalar*(a*b*uetp);
T(2).d = -scalar*(c*uetp);
T(3).d = -scalar*(a*uerp);
T(4).d = -scalar*(b*uerp);
T(5).d = -scalar*(d*uetp);
