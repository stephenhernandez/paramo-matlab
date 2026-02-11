function [R,rDiff,X,Y,dobs,B,d1,d2,d1d2,synth,sizeX,S] = sierraNegraTremorAmplitudeLocation(S)
f = 1.75; %center frequency
Q = 40; %attenuation quality factor
beta = 2; %11/10; %seismic shear wave velocity
B = pi*f/(Q*beta);

%%
hInc = 0.002;
[~,~,minLat,maxLat,minLon,maxLon] = get_region_dimensions('sierra_negra');
minLon = minLon+0.3;
maxLon = maxLon-0.3;
minLat = minLat + 0.3;
maxLat = maxLat + 0.3;

%%
[X,Y] = meshgrid(minLon:hInc:maxLon,minLat:hInc:maxLat);
sizeX = size(X);
X = X(:);
Y = Y(:);
lX = length(X);
refEllipse = referenceEllipsoid('wgs84');
depth = 0;

%%
filtN = 21;
S = medfiltWaveforms(S,filtN);
S = filterWaveforms(S,-inf,1/(S(1).delta*filtN),1,[],true);

S = syncWaveforms(S);
lS = length(S);
stla = NaN(lS,1);
stlo = stla;
for i = 1:lS
    kstnm = S(i).kstnm;
    [sensorCode,digitizerCode] = getResponseCodes(kstnm);
    digitizerGain = readDigitizerGain(digitizerCode);
    [~,~,sensitivity,~] = readSensorPolesZeros(sensorCode);
    constant_ = 1./(sensitivity*digitizerGain);
    S(i) = scaleWaveforms(S(i),1e6*constant_);
    [stla_,stlo_] = metaDataFromStationList(kstnm);
    stla(i) = stla_;
    stlo(i) = stlo_;
    %     if strcmp(kstnm,"PVIL")
    %         S(i) = scaleWaveforms(S(i),true);
    %     elseif strcmp(kstnm,"CEAZ")
    %         S(i) = scaleWaveforms(S(i),true);
    %     end
    if strcmp(kstnm,"ALCE")
        S(i) = scaleWaveforms(S(i),1/8);
    end
end

%% get ratios
Ndiff = 0.5*lS*(lS-1);
npts = S(1).npts;
R = repmat(S(1),Ndiff,1);   % observed ratio
rDiff = NaN(lX,npts,Ndiff); % synthetic - observed
n = 0;
for i = 1:lS-1
    d1 = distance(Y,X,stla(i),stlo(i),refEllipse)*1e-3;
    d1 = sqrt(d1.^2 + depth^2);
    for j = i+1:lS
        n = n+1;
        R(n) = S(i);
        R(n).d = S(i).d./S(j).d;
        %R(n).d = R(n).d./R(n).d;
        %R(n) = medfiltWaveforms(R(n),filtN);
        %R(n) = filterWaveforms(R(n),-inf,1/(R(n).delta*filtN),2,[],true);
        d2 = distance(Y,X,stla(j),stlo(j),refEllipse)*1e-3;
        d2 = sqrt(d2.^2 + depth^2);
        d1d2 = d1-d2;
        %d2 = d2./d1; %
        d2 = sqrt(d2./d1);
        dobs = R(n).d;
        dobs = dobs(:)';
        synth = d2.*exp(-B.*d1d2);
        
        for k = 1:npts
            r_ = (synth - dobs(k)'); %./synth;
            r_ = abs(r_);
            r_ = r_/sum(r_);
            rDiff(:,k,n) = r_;
        end
    end
end