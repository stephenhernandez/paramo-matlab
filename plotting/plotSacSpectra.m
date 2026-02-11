function ax = plotSacSpectra(S,nfft,diffFlag,lfc,hfc)
if nargin < 5; hfc = -inf; end
if nargin < 4; lfc = -inf; end
if nargin < 3; diffFlag = false; end
if nargin < 2; nfft = 2048; end

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
    zeroPhaseFlag = true;
    data = detrend(diff(data));
    data = zpkFilter(data,lfc,hfc,Fs,npoles,zeroPhaseFlag);
    data = cumsum(data);
else
    disp('not filtering');
end
ld = length(data);
tdum = (0:ld-1)/Fs;

figure('units','normalized','outerposition',[0 0 1 1]);
ax = gobjects(Nchan,2);
for i = 1:Nchan
    disp(i)
    ax(i,1) = subplot(Nchan,2,2*(i-1)+1);
    plot(ax(i,1),tdum,data(:,i),'linewidth',2);
    title(strcat(S(i).kstnm,',',S(i).kcmpnm));
    ax(i,1).XTickLabelRotation = 15;
    axis(ax(i,1),'tight');
    ax(i,1).YLabel.String = 'Amp. [Counts]';
    ax(i,1).XLabel.String = ['Seconds since: ',datestr(S(1).ref)];
    
    ax(i,2) = subplot(Nchan,2,2*i);
    [pxx,fxx] = pmtm(data(:,i),3,nfft,Fs);
    semilogx(ax(i,2),fxx,pxx,'linewidth',2);
    %plot(ax(i,2),fxx,pxx);
    axis(ax(i,2),'tight');
    fi_lf = sum(pxx(fxx >= 0.2 & fxx < 0.5));
    fi_hf = sum(pxx(fxx >= 0.5 & fxx < 5));
    fi = fi_lf/fi_hf;
    ax(i,2).YLabel.String = 'Power';
    %ax(i,2).XLabel.String = 'Period [Sec.]';
    ax(i,2).XLabel.String = 'Frequency [Hz.]';
    ax(i,2).XLim = [0.1 25];
    %legend(['Low Band/High Band: ',num2str(fi)],'location','northwest');
end
linkaxes(ax(:,1),'x');
linkaxes(ax(:,2),'x');

