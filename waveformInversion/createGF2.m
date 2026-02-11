%function g = createGF2(az,dist,filtflag,Hd,stf,starttime,duration,depth,flag3,synthflag,Fs_obs,stnm,unitFlag)
function g = createGF2(baz,az,dist,filtflag,lfc,hfc,npoles,stf,starttime,duration,depth,flag3,synthflag,Fs_obs,stnm,unitFlag)

%dist should be in km
if mod(depth,1) == 0
    depth = floor(depth) + 0.5;
end

%
dist = round(dist)*10;
gstr = num2str(dist);
if dist < 10
    gstr = strcat('000',gstr);
elseif dist < 100
    gstr = strcat('00',gstr);
elseif dist < 1000
    gstr= strcat('0',gstr);
end

%
if synthflag
    %prefix = '/auto/proj/stephenh/amplitudes/nicoyaGF/';
    %prefix = '/Volumes/HERNARADO_Toshiba_2/gf/'; %for external HD
    %prefix = '~/gf/nicoyaGF/'; %for external HD
    prefix = '~/gf/ecuador/agudelo1Hz_triangleSTF_L2_Disp/';
    %prefix = '~/gf/ecuador/vallee/'; %for external HD (VR reduction not as high as CR)
    %prefix = '~/gf/costa_rica/agudelo/'; %for external HD
    %prefix = '~/gf/ecuador/agudelo/'; %for external HD
    %prefix = '~/gf/ecuador/ASW2/'; %for external HD
else
    prefix = '~/gf/GFS_0.1d/';
end

dpstr = num2str(depth);
if depth < 10
    dpstr = strcat('0',dpstr);
    % elseif dist < 100
    %     dpstr= strcat('0',dpstr);
end
dpstr = strcat('H0',dpstr,'/');

base = strcat('GF.',gstr,'.SY.');
az = deg2rad(az);

%%
try
    [uzrr,h] = readsacfile(strcat(prefix,dpstr,'RR/',base,'LHZ.SAC'));  %figure; plot(uzrr,'k');
    uztt = readsacfile(strcat(prefix,dpstr,'TT/',base,'LHZ.SAC'));
    uzpp = readsacfile(strcat(prefix,dpstr,'PP/',base,'LHZ.SAC'));
    uzrt = readsacfile(strcat(prefix,dpstr,'RT/',base,'LHZ.SAC'));
    usrr = readsacfile(strcat(prefix,dpstr,'RR/',base,'LHL.SAC'));
    ustt = readsacfile(strcat(prefix,dpstr,'TT/',base,'LHL.SAC'));
    uspp = readsacfile(strcat(prefix,dpstr,'PP/',base,'LHL.SAC'));
    usrt = readsacfile(strcat(prefix,dpstr,'RT/',base,'LHL.SAC'));
    uerp = readsacfile(strcat(prefix,dpstr,'RP/',base,'LHT.SAC'));
    uetp = readsacfile(strcat(prefix,dpstr,'TP/',base,'LHT.SAC'));
catch ME
    fprintf('%s\n',strcat(prefix,dpstr,'RR/',base,'LHZ.SAC'));
    rethrow(ME);
end

%%
dt = h(1);
Fs_syn = round(1/dt);
Hd = zpkOperator(lfc,hfc,Fs_syn,npoles);

convopt = 'full';
a = sin(az); 
b = cos(az);
c = sin(2*az); 
d = cos(2*az);

a2 = a*a;
b2 = b*b;
a2b2 = a2 - b2;

zx1 = uzrr - a2*uztt - b2*uzpp;
zx2 = d*uztt + a2b2*uzpp;
zx3 = b*uzrt;
zx4 = a*uzrt;
zx5 = c*(uztt-uzpp);

rx1 = (usrr - a2*ustt - b2*uspp);
rx2 = (d*ustt + a2b2*uspp);
rx3 = b*usrt;
rx4 = a*usrt;
rx5 = c*(ustt-uspp);

tx1 = a*b*uetp;
tx2 = c*uetp;
tx3 = a*uerp;
tx4 = b*uerp;
tx5 = d*uetp;

Z1 = [zx1 zx2 zx3 zx4 zx5];
R1 = [rx1 rx2 rx3 rx4 rx5];
T1 = [tx1 tx2 tx3 tx4 tx5];

if strcmp(unitFlag, 'acc')
    % convert to acceleration
    Z1 = diff(diff(Z1)/dt)/dt;
    R1 = diff(diff(R1)/dt)/dt;
    T1 = diff(diff(T1)/dt)/dt;
elseif strcmp(unitFlag, 'vel')
    % convert to velocity
    Z1 = diff(Z1)/dt;
    R1 = diff(R1)/dt;
    T1 = diff(T1)/dt;
end

if Fs_syn ~= Fs_obs
    disp('synthetics and observables have dissimilar sample rate. quitting')
    return
end

if filtflag
    %convolve
    Z1 = convn(Z1,stf,convopt);
    R1 = convn(R1,stf,convopt);
    T1 = convn(T1,stf,convopt);

    %filter (must come before any cutting)
    Z1 = filter(Hd,Z1);
    R1 = filter(Hd,R1);
    T1 = filter(Hd,T1);

    t = (0:length(Z1(:,1))-1)*dt;
    t = shiftdim(t);
    tI = t >= starttime;
    
    Z1t = zeros(sum(tI),5); 
    R1t = Z1t; 
    T1t = Z1t;

    %
    for j = 1:5
        %cut the data
        Ztmp = Z1(:,j);
        Z1t(:,j) = Ztmp(tI);
        Rtmp = R1(:,j);
        R1t(:,j) = Rtmp(tI);
        Ttmp = T1(:,j);
        T1t(:,j) = Ttmp(tI);
    end    
else
    
    t = (0:length(Z1(:,1))-1)*dt;
    t = shiftdim(t);
    tI = t >= starttime;% & t < endtime;
    Z1t = zeros(sum(tI),5); R1t = Z1t; T1t = Z1t;
    for j = 1:5
        %cut the data
        Ztmp = Z1(:,j);
        Z1t(:,j) = Ztmp(tI);
        Rtmp = R1(:,j);
        R1t(:,j) = Rtmp(tI);
        Ttmp = T1(:,j);
        T1t(:,j) = Ttmp(tI);
    end
end

%%
rotang = baz;
if baz <= 180
    rotang = baz + 180;
elseif baz > 180 && baz < 360
    rotang = baz - 180;
end

for ii = 1:5
    rtemp = R1t(:,ii);
    ttemp = T1t(:,ii);
    [ee,nn] = rotate2d(ttemp,rtemp,rotang);
    R1t(:,ii) = nn;
    T1t(:,ii) = ee;
end

if flag3
    %g = -[Z1t(1:duration,:); R1t(1:duration,:); T1t(1:duration,:)];
    %g = -[T1t(1:duration,:); R1t(1:duration,:); Z1t(1:duration,:)];
    %g = -[zeros(duration,5); R1t(1:duration,:); Z1t(1:duration,:)];
    if strcmp(stnm,'ACHA1')
        %disp(stnm)
        g = [zeros(duration,5); R1t(1:duration,:); Z1t(1:duration,:)];
    elseif strcmp(stnm,'POPE1')
        g = [zeros(duration,5); R1t(1:duration,:); Z1t(1:duration,:)];
    else
        g = [T1t(1:duration,:); -R1t(1:duration,:); -Z1t(1:duration,:)];
    end
else
    %     %just the horizontals
    %     if strcmp(stnm,'POPE1')
    %         disp(stnm)
    %         g = [zeros(duration,5); R1t(1:duration,:)];
    %     elseif strcmp(stnm,'ACHA1')
    %         g = [zeros(duration,5); R1t(1:duration,:)];
    %     else
    %         %g = [T1t(1:duration,:); R1t(1:duration,:)];
    %g = [R1t(1:duration,:); Z1t(1:duration,:)];
    %     end
    
    %g = T1t(1:duration,:);
    %g = R1t(1:duration,:);
    g = Z1t(1:duration,:);
end
