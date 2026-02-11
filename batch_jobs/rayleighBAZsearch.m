function [estBAZ,corrcoefs,angs,corrCoeffEstBAZ,loveTT,rayTT,radial,transverse] = rayleighBAZsearch(S,lfc,hfc,plotFlag,origBAZ)
if nargin < 2; lfc = -inf; end
if nargin < 3; hfc = -inf; end
if nargin < 4; plotFlag = false; end
if nargin < 5; origBAZ = NaN; end

%npoles = 4;
S = syncWaveforms(S);
%S = transferWaveforms(S,lfc,hfc,npoles,1,'vel');
%S = scaleWaveforms(S,1e3);

%%
%Fs = round(1/S(1).delta);
delta = median(pull(S,'delta'));
%% filter the data
if any(isfinite([lfc hfc]))
    disp('filtering data')
    disp([lfc hfc])
    tw = 0.005;
    S = filterWaveforms(taperWaveforms(detrendWaveforms(S),2./lfc./delta),lfc,hfc);
else

end
data = pull(S);
data = detrend(data);

npts=size(data,1);
e = data(1:npts,3);
n = data(1:npts,2);
z = data(1:npts,1);

%zeroPhaseFlag = false;
ef = e; %taper(detrend(zpkFilter(e,lfc,hfc,Fs,npoles,zeroPhaseFlag)),tw);
nf = n; %taper(detrend(zpkFilter(n,lfc,hfc,Fs,npoles,zeroPhaseFlag)),tw);
zf = z;



zh=hilbert(zf);
zfh=imag(zh);
normzfh = norm(zfh);
zfh = zfh./normzfh;

inc = 1;
angs = 0:inc:359;
Nangs = length(angs);
corrcoefs = zeros(Nangs,1);
for i = 1:Nangs
    theta = angs(i);
    [erot,nrot] = rotate2d(ef,nf,theta);
    %nrot = rotate2d(ef,nf,theta);
    dum_ = dot(nrot,zfh)./norm(nrot); %corrcoef([nrot,zfh]);
    dum2_ = dot(erot,zfh)./norm(erot); %corrcoef([nrot,zfh]);
    corrcoefs(i) = dum_(1); %,2);
    cc2(i) = dum2_(1);
end
[corrCoeffEstBAZ,cI] = max(corrcoefs);
estBAZ = angs(cI);

t = getTimeVec(S(1));
[erot,nrot] = rotate2d(ef,nf,estBAZ);
%t = seconds(t-min(t));
erotEnv = abs(hilbert(erot));
nrotEnv = abs(hilbert(nrot));
[~,mI1] = max(erotEnv);
[~,mI2] = max(nrotEnv);
loveTT = t(mI1);
rayTT = t(mI2);

% radial = peak2peak(nrot);
% transverse = peak2peak(erot);

radial = sqrt(detrend(nrot).^2 + detrend(zf).^2);
transverse = detrend(erot);

if plotFlag
    figure('units','normalized','outerposition',[0 0 1 1]);
    ha(1) = subplot(411);
    plot(t,erot); %hold on; plot(t,normzfh.*zfh);
    ylabel('transverse'); axis tight;
    title(sprintf('%s.%s.%s',S(1).knetwk,S(1).khole,S(1).kstnm));
    ha(2) = subplot(412);
    plot(t,nrot); hold on; plot(t,normzfh.*zfh); legend('R','h{Z}');
    ylabel('radial'); axis tight;
    ha(3) = subplot(413);
    plot(t,zf); ylabel('vertical'); axis tight;
    linkaxes(ha,'x');
    zoom on;

    subplot(414);
    if isfinite(origBAZ)
        [maxCorr,mI] = max(corrcoefs);
        plot(angs-origBAZ,corrcoefs,'o');
        title(strcat('Orig BAZ:',num2str(origBAZ),',',num2str(maxCorr),',',...
            num2str(round(angs(mI)-origBAZ))));
        xlim([0 360]);
    else
        plot(angs,corrcoefs,'o');
        hold on;
        plot(angs,cc2,'s');
        title(['Estimated BAZ: ',num2str(estBAZ),', Best CC: ',num2str(corrCoeffEstBAZ)]);
        xlim([0 360]);
        plot(angs(cI),corrCoeffEstBAZ,'p','markerfacecolor','w','markeredgecolor','k','markersize',15);
        legend('R','T','best solution');
    end
end
