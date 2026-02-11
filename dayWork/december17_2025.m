%17 DEC 2025
%testing a rotation algorithm to maximize SNR on a small aperture array, to
%then determine doa and app. vel.
clear; close all;
kstnms = ["CTX1";"CTX2";"CTX3";"CTX4";"CTX5";"CTX6";"CTX7"]; transferFlag = false;
elementLats = [-0.65787; -0.65682; -0.65602; -0.65567; -0.65543; -0.65421; -0.65452];
elementLons = [-78.43968; -78.43797; -78.44049; -78.43839; -78.43647; -78.43969; -78.43798];
elementElevs = [4643; 4625; 4581; 4587; 4551; 4542; 4534];
rawdatadir = "~/data/nodes/CotopaxiNodeData2023";
lfc = 0.6; hfc = 1.2;
newFs = 2000;
for i = 1:length(kstnms)
    kstnm_ = kstnms(i);
    C = loadWaveforms(datetime(2023,05,12),1,kstnm_,["DPZ";"DPN";"DPE"],"EC","",true,true,rawdatadir);
    Ccut = detrendWaveforms(cutWaveforms(C,dateshift(C(1).ref,'start','day')+hours(13)+minutes(04)+seconds(00),0,minutes(1)));
    
    Ccut = syncWaveforms(Ccut,0,1,1);
    Ccut = resampleWaveforms(Ccut,newFs);
    Cf = filterWaveforms(taperWaveforms(detrendWaveforms(Ccut),newFs*10),lfc,hfc);
    %Cf(1) = dealHeader(Cf(1),imag(hilbert(pull(Cf(1)))));
    d = pull(Cf);
    N = size(d,1);
    [u1,s1,v1] = svd(d(1:N,:),"econ");
    az_ = 90-atan2d(v1(2,1),v1(3,1))
    Cf = rotateWaveforms(Cf,az_,3,2,true);
    d = pull(Cf);
    N = size(d,1);
    [u,s,v] = svd(d(1:N,:),"econ");
    inc_ = 90-atan2d(v(1,1),v(2,1))
    Cf2 = rotateWaveforms(Cf,inc_,2,1,true);
    d2 = pull(Cf2);
    N = size(d2,1);
    [u2,s2,v2] = svd(d2(1:N,:),"econ");
    Csave(i,1) = Cf(1);
    Csave(i) = dealHeader(Csave(i),u2(:,1));
end
%
d = pull(Csave);
N = size(d,1);
[u,s,v] = svd(d(1:N,:),"econ"); svOrig = s*v';
close all; [~,ax1] = plotWaveforms(Csave);
N = 4; 
R = u(:,1:N)*s(1:N,1:N)*v(:,1:N)';
Cf1 = Csave; Cf2 = Csave;
for i = 1:7
    v_ = svOrig(:,i).^2;
    [~,sI] = sort(v_,"descend");
    disp(sI')
    Cf1(i) = dealHeader(Cf1(i),u(:,sI(1:N))*svOrig(sI(1:N),i));
    Cf2(i) = dealHeader(Cf2(i),R(:,i));
end
plotWaveforms(ax1,Cf1);
%plotWaveforms(ax1,Cf2);
res = Csave; d = pull(Csave); dpred = pull(Cf1); res_ = d-dpred; for i = 1:7; res(i) = dealHeader(res(i),res_(:,i)); end; plotWaveforms(res);
%
Cf1 = normalizeWaveforms(Cf1); %denoised
[u,s,v] = svd(pull(Cf1),"econ");
s
v
s*v'
sum((s*v').^2)
rssq(pull(Cf1))
cumsum((s*v').^2)
svOrig.^2
dpred = pull(fdWhiten(Cf1,lfc,hfc));
[maxccp,plags] = doccFreqCircShiftPolarities(dpred,true,[],newFs/2);
[shifted_data,maxccp_,G,plags_,raw_shifts,meancc,medcc] = apply_vdcc(dpred,[maxccp,plags]);
(max(raw_shifts)-raw_shifts)/newFs
linkedPlot([],shifted_data,"-");
refEllipse = referenceEllipsoid('wgs84');
lS = length(elementLats);
[d_,az_] = distance(mean(elementLats),mean(elementLons),...
    elementLats,elementLons,refEllipse);
Gorig = [d_.*cosd(90-az_),d_.*sind(90-az_),elementElevs - mean(elementElevs)];
Gorig = -[getDD(Gorig(:,1)) getDD(Gorig(:,2))];% getDD(Gorig(:,3))];
Gcolumns = size(Gorig,2);
Ginv = pinv(Gorig);
[doa_,av_] = doa_av_inversion(-getDD((max(raw_shifts)-raw_shifts))/newFs,Ginv)
