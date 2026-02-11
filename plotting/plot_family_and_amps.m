function tmp_stack = plot_family_and_amps(indiv_events,t,amps,indices,ampFact,Fs,lfc,hfc,zeroPhaseFlag,tbounds,abounds)
if nargin < 5; ampFact = 5; end
if nargin < 6; Fs = 100; end
if nargin < 7; lfc = -inf; end
if nargin < 8; hfc = -inf; end
if nargin < 9; zeroPhaseFlag = false; end
if nargin < 12; hFlag = false; end
% if nargin < 10; tbounds = []; end %datetime(2016,04,01),datetime(2016,04,60)]; end
% if nargin < 11; abounds = [10^floor(log10(min(amps))) 10^ceil(log10(max(amps)))]; end %10 1.5e4]; end

npoles = 4;
if any(isfinite([lfc hfc]))
    disp('filtering data')
    indiv_events = zpkFilt(indiv_events,lfc,hfc,Fs,npoles,zeroPhaseFlag);
else
    disp('no filtering requested');
end
indiv_tmp = indiv_events(:,indices);

%%
[winlen,nEvents] = size(indiv_tmp);
indiv_tmp = -normalizeWaveforms(indiv_tmp); %flip here to offset flipping done by ax.YDir = 'reverse';
if nEvents >= 10
    %tmp_stack = nanmedian(indiv_tmp,2); %
    tmp_stack = pws(normalizeWaveforms(detrend(indiv_tmp))); %phase-weighted stack
else
    tmp_stack = mean(indiv_tmp,2,"omitnan"); %
    %tmp_stack = pws(normalizeWaveforms(detrend(indiv_tmp))); %phase-weighted stack
end

tmp_stack = tmp_stack/norm(tmp_stack);
tdum = (0:winlen-1)/Fs;
li = length(indices);

for i = 1:li
    indiv_tmp(:,i) = ampFact*indiv_tmp(:,i)+i-1;
end
indiv_tmp = [indiv_tmp ampFact*tmp_stack+li];

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(4,3,[2 3 5 6 8 9]);
ax = gca;
plot(ax,t,amps,'o','linewidth',2);
grid on;
ylabel('amplitude');
if nargin >= 10
    set_xlim(ax,tbounds)
    set_ylim(ax,abounds);
end
ax.YScale = 'log';
%xtick_ = dateshift(linspace(tbounds(1),tbounds(2),12),'start','day','nearest');
%ax.XTick = xtick_;

subplot(4,3,[1 4 7]);
ax = gca;
plot(ax,tdum,indiv_tmp,'linewidth',2);
hold(ax,'on');
plot(ax,tdum,indiv_tmp(:,end),'k','linewidth',3);
ax.YDir = 'reverse';
axis tight
grid on;
xlabel('Seconds');
ylabel('Event Number')

subplot(4,3,[10 11 12]);
plot(tdum,-normalizeWaveforms(indiv_tmp,1,1),'-','color',[0.5 0.5 0.5],'linewidth',2);
hold on;
plot(tdum,-normalizeWaveforms(tmp_stack,1,1),'k','linewidth',2);
zoom on;

if hFlag
    plot(tdum,abs(hilbert(-normalizeWaveforms(tmp_stack,1,1))),'r','linewidth',3);
end
xlabel('time [sec.]')
ylabel('Normalized Amplitude');
axis tight;
tmp_stack = -normalizeWaveforms(detrend(tmp_stack),1);

function set_xlim(ax,vlines)
lvlines = length(vlines);
if lvlines
    ax.XLim = vlines;
else
    return
end

function set_ylim(ax,ybounds)
lvlines = length(ybounds);
if lvlines
    ax.YLim = ybounds;
else
    return
end