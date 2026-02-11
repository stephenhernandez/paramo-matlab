function ax = plotSacSpectraSpectrograms(S,window,noverlap,nfft,diffFlag,lfc,hfc,satLev,floorLev)
if nargin < 9; floorLev = 0; end
if nargin < 8; satLev = 50; end
if nargin < 7; hfc = -inf; end
if nargin < 6; lfc = -inf; end
if nargin < 5; diffFlag = false; end
if nargin < 4; nfft = 2048; end
if nargin < 3; noverlap = 512; end
if nargin < 2; window = 1024; end

if diffFlag
    S = differentiateSacData(S);
end

Nchan = length(S);
if Nchan > 1
    S = synchSacData(S);
end
Fs = 1/S(1).delta;
data = pull(S);

if any(isfinite([lfc hfc]))
    disp('filtering data')
    disp([lfc hfc]);
    npoles = 4;
    zeroPhaseFlag = false;
    data = detrend(diff(data));
    data = zpkFilter(data,lfc,hfc,Fs,npoles,zeroPhaseFlag);
    data = cumsum(data);
end
ld = length(data);
tdum = (0:ld-1)/Fs;

figure('units','normalized','outerposition',[0 0 1 1]);
ax = gobjects(2,Nchan);
for i = 1:Nchan
    disp(i);
    ax(1,i) = subplot(2,Nchan,i);
    plot(tdum,data(:,i),'linewidth',2);
    title(strcat(S(i).kstnm,',',S(i).kcmpnm));
    ax(1,i) = gca;
    ax(1,i).XTickLabelRotation = 15;
    axis tight;
    hc = colorbar;
    hc.Visible = 'off';
    
    ax(2,i) = subplot(2,Nchan,i+Nchan);
    [~,ff,tt,pp] = spectrogram(S(i).d,parzenwin(window),noverlap,nfft,Fs);
    %[~,ff,tt,pp] = spectrogram(data(:,i),window,noverlap,nfft,Fs);
    fI = ff > 0;
    dataSpec = 10*log10(abs(pp));
    imagesc(tt,ff(fI),dataSpec(fI,:)); 
    hold on;
    axis xy;
    colorbar;
    caxis([floorLev satLev]);
    ylabel('Frequency [Hz.]');
    axis tight;
    ax(2,i) = gca;
    ax(2,i).XTickLabelRotation = 15;
    ax(2,i).XLabel.String = ['Seconds since: ',datestr(S(1).ref)];
    %ax(2,i).YScale = 'log';
    ylim([0 25]);
    cmap = jet(512);
    colormap(cmap);
end
linkaxes(ax,'x');
linkaxes(ax(2,:),'xy');