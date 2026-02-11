function S = getSyntheticTRZ(az,dist,depth,strike,dip,rake,newFs,bazFlag)
if nargin < 8
    bazFlag = true;
end

if nargin < 7
    newFs = 2;
end

if mod(depth,1) == 0
    depth = depth + 0.5;
end

dist = round(dist)*10;
gstr = num2str(dist);
if dist < 100
    gstr = strcat('00',gstr);
elseif dist < 1000
    gstr= strcat('0',gstr);
end

prefix = '~/gf/nicoyaGF/';

dpstr = num2str(depth);
if depth < 10
    dpstr = strcat('0',dpstr);
end
dpstr = strcat('H0',dpstr,'/');

base = strcat('GF.',gstr,'.SY.');
az = deg2rad(az);

uzrrFile = strcat(prefix,dpstr,'RR/',base,'LHZ.SAC');
disp(uzrrFile);

%%
S = sac2struct(uzrrFile);
S(2) = S(1);
S(3) = S(1);

%%
uzrr = S(1).d;
uztt = readsacfile(strcat(prefix,dpstr,'TT/',base,'LHZ.SAC'));
uzpp = readsacfile(strcat(prefix,dpstr,'PP/',base,'LHZ.SAC'));
uzrt = readsacfile(strcat(prefix,dpstr,'RT/',base,'LHZ.SAC'));
usrr = readsacfile(strcat(prefix,dpstr,'RR/',base,'LHL.SAC'));
ustt = readsacfile(strcat(prefix,dpstr,'TT/',base,'LHL.SAC'));
uspp = readsacfile(strcat(prefix,dpstr,'PP/',base,'LHL.SAC'));
usrt = readsacfile(strcat(prefix,dpstr,'RT/',base,'LHL.SAC'));
uerp = readsacfile(strcat(prefix,dpstr,'RP/',base,'LHT.SAC'));
uetp = readsacfile(strcat(prefix,dpstr,'TP/',base,'LHT.SAC'));

%%
a = sin(az); 
b = cos(az);
c = sin(2*az); 
d = cos(2*az);

a2 = a*a;
b2 = b*b;
a2b2 = a2 - b2;

%%
zx1 = uzrr - a2*uztt - b2*uzpp;
zx2 = d*uztt + a2b2*uzpp;
zx3 = b*uzrt;
zx4 = -a*uzrt;
zx5 = -c*(uztt-uzpp);

rx1 = usrr - a2*ustt - b2*uspp;
rx2 = d*ustt + a2b2*uspp;
rx3 = b*usrt;
rx4 = -a*usrt;
rx5 = -c*(ustt-uspp);

tx1 = -a*b*uetp;
tx2 = -c*uetp;
tx3 = -a*uerp;
tx4 = -b*uerp;
tx5 = -d*uetp;

%%
Z1 = [zx1 zx2 zx3 zx4 zx5];
R1 = [rx1 rx2 rx3 rx4 rx5];
T1 = [tx1 tx2 tx3 tx4 tx5];

%%
sdr = [strike,dip,rake];
mt = sdr2mt(sdr)';
mt(3) = [];

%%
Z1 = Z1*mt;
R1 = R1*mt;
T1 = T1*mt;

%%
S(1).kcmpnm = "T";
S(2).kcmpnm = "R";
S(3).kcmpnm = "Z";

%%
if bazFlag
    baz = rad2deg(az) + 180;
    if baz >= 360
        baz = baz - 360;
    end
    
    if baz <= 180
        rotang = baz + 180;
    else
        rotang = baz - 180;
    end
    rotang = -rotang;
    
    disp(rotang)
    [T1,R1] = rotate2d(T1,R1,rotang);
    S(1).kcmpnm = "E";
    S(2).kcmpnm = "N";
end

%%
S(1).d = T1;
S(2).d = R1;
S(3).d = Z1;

%%
S = resampleWaveforms(S,newFs);
