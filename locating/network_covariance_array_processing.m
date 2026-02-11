function [DOA,AV,DT] = network_covariance_array_processing(CCpairs,lfc,hfc,...
    elementLats,elementLons,elementElevs,newFs,npoles)
[~,ntime,npairs] = size(CCpairs);
upsampleFactor = 100;
Fs2 = newFs*upsampleFactor;

%%
refEllipse = referenceEllipsoid('wgs84');
lS = length(elementLats);

[d_,az_] = distance(mean(elementLats),mean(elementLons),...
    elementLats,elementLons,refEllipse);
Gorig = [d_.*cosd(90-az_),d_.*sind(90-az_),elementElevs - mean(elementElevs)];
Gorig = -[getDD(Gorig(:,1)) getDD(Gorig(:,2))];% getDD(Gorig(:,3))];
Gcolumns = size(Gorig,2);

Ginv = pinv(Gorig);
%%
lStations = (-1 + sqrt(1 + 4*(2*npairs)))/2;
if lS ~= lStations
    %complain hard
    fprintf('number of stations dont match, fatal error\n');
    return;
end

%%
nXpairs = lStations*(lStations-1)*0.5;
% autocorr_pairs = 1+cumsum([0;(lStations:-1:1)']);
% autocorr_pairs = autocorr_pairs(1:lStations);
autocorr_pairs = (nXpairs+1:npairs)';
DOA = NaN(ntime,1);
AV = DOA;
%DT_DIFF = DOA;
if Gcolumns > 2
    ELEVANGLE = DOA;
else
    ELEVANGLE = [];
end

DT = NaN(nXpairs,ntime);
%Gshift = Gvdcc(lS);

%Fix this section, do it in frequency domain, make it zero phase by
%default...
Hd = zpkOperator(lfc,hfc,newFs,npoles);
df = ifft(CCpairs,[],1,'symmetric'); %(:,i,:);
clear CCpairs;
df = fftshift(df,1);
df(:,:,autocorr_pairs) = [];

df = filter(Hd,df);
df = flipud(df);
df = filter(Hd,df);
df = flipud(df);
df_ = resample(squeeze(df(:,1,1)),upsampleFactor,1);
m = size(df_,1)/2;
lags = (-m:m-1)';

for i = 1:ntime
    fprintf("%d\n",i);
    df_ = squeeze(df(:,i,:));
    df_ = resample(df_,upsampleFactor,1);

    [~,maxI] = max(df_,[],1);
    dT = lags(maxI);
    dT = dT / Fs2;
    DT(:,i) = dT;

    if Gcolumns > 2
        [doa_,av_,elev_angle_] = doa_av_inversion(dT,Ginv);
        ELEVANGLE(i) = elev_angle_;
    else
        [doa_,av_] = doa_av_inversion(dT,Ginv);
    end
    DOA(i) = doa_;
    AV(i) = av_;
end