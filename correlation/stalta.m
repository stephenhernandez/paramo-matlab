function [locs,snr,staOlta,sosSTA] = stalta(S,sta,lta,mph,hFlag,plotFlag,envFiltFlag,hfc,verboseFlag)

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
%[locs,snr,staOlta,sosSTA] = stalta(S,sta,lta,mph,hFlag,plotFlag,envFiltFlag,hfc)
if nargin < 2; sta = 5; end
if nargin < 3; lta = 10; end
if nargin < 4; mph = 3; end
if nargin < 5; hFlag = false; end
if nargin < 6; plotFlag = false; end
if nargin < 7; envFiltFlag = false; end
if nargin < 8; hfc = 1/sta; end
if nargin < 9; verboseFlag = false; end

%%
S = S(1);
gapFlag = S.gapFlag;
if gapFlag
    gstops = sum(S.gapInfo,2) - 1;
else
    gstops = [];
end

%%
delta = S.delta;
if delta < 1
    Fs = round(1/delta);
else
    Fs = 1/delta;
end
n_sta = round(sta*Fs);
n_lta = round(lta*Fs);

%%
if verboseFlag
    fprintf('Using STA = %d points\n',n_sta);
    fprintf('Using LTA = %d points\n',n_lta);
end

%%
d = S.d;
if hFlag
    %disp('Getting envelope')
    d2 = abs(hilbert(d));
else
    d2 = abs(d).^2;
end
%zeroI = d2 <= 1;

%%
boxSTA = ones(n_sta,1)/n_sta;
if n_sta > 1
    sosSTA = flipud(fftfilt(boxSTA,flipud(d2)));
elseif n_sta == 1
    sosSTA = d2;
end

%%
if n_lta < 10
    sosLTA = medfiltSH(d2,n_lta);
else
    boxLTA = ones(n_lta,1)/n_lta;
    sosLTA = fftfilt(boxLTA,d2);
end

%%
staOlta = sqrt(abs(sosSTA))./sqrt(abs(sosLTA));
staOlta(1:n_lta) = 1;
zeroI = staOlta <= 1;
staOlta(zeroI) = 0;

%%
if envFiltFlag
    disp('smoothing the envelope');
    staOlta = zpkFilter(staOlta,-inf,hfc,Fs,1,1);
end

%%
[snr,locs]= findpeaks(staOlta,'MinPeakHeight',mph,'MinPeakDistance',1.*n_sta);
locsI = locs >= n_lta & isfinite(snr);
locs = locs(locsI);
snr = snr(locsI);

%%
nDetections = sum(locsI);
if ~nDetections
    if plotFlag
        fprintf('no events found, not plotting anything\n');
    else
        if verboseFlag
            fprintf('no events found\n');
        end
    end
    return;
end

if ~isempty(gstops)
    trueI = true(length(locs),1);
    thresh = n_sta;
    for i = 1:length(locs)
        if min(abs(locs(i)-gstops)) <= thresh
            trueI(i) = false;
        end
    end
    locs = locs(trueI);
    snr = snr(trueI);
end

if plotFlag
    t = getTimeVec(S);
    lt = length(t);
    l = min([lt length(staOlta)]);

    figure('units','normalized','outerposition',[0 0 1 1]);
    ha(1) = subplot(311);
    plot(t(1:l),sosSTA); hold on;
    plot(t(1:l),sosLTA);
    plot(t(1:l),d2(1:l));
    legend('sta','lta','envelope','location','northwest');
    set(gca,'YScale','log');
    ha(2) = subplot(312);

    plot(t(1:l),staOlta(1:l)); hold on;

    pp = plot(t(locs),snr,'o');
%    pp.LineWidth = 1;
    ylabel('snr');
    set(gca,'YScale','log');
    grid on;

    ha(3) = subplot(313);
    plot(t(1:l),d(1:l));
    zoom on;
    linkaxes(ha,'x');
end
