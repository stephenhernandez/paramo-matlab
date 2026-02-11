function tmp_stack = plot_family(indiv_events,indices,ampFact,Fs,lfc,hfc,zeroPhaseFlag,stackFlag,shift)
if nargin < 4; Fs = 100; end
if nargin < 5; lfc = -inf; end
if nargin < 6; hfc = -inf; end
if nargin < 7; zeroPhaseFlag = false; end
if nargin < 8; stackFlag = true; end
if nargin < 9; shift = 0; end

if any(isfinite([lfc hfc]))
    npoles = 4;
    fprintf("filtering data\n");
    tw = ceil(2*Fs/lfc);
    indiv_events = diff(indiv_events);
    indiv_events = taper(indiv_events,tw);
    indiv_events = zpkFilter(indiv_events,lfc,hfc,Fs,npoles,zeroPhaseFlag);
    indiv_events = cumsum(indiv_events);
else
    fprintf("no filtering requested\n");
end
indiv_tmp = indiv_events(:,indices);

%%
indiv_tmp = -normalizeWaveforms(indiv_tmp); %flip here to offset flipping done by ax.YDir = 'reverse';
%[U,~,V] = svd(normalizeWaveforms(detrend(indiv_tmp)),"econ");
%U = U(:,1);
%%mean_pol = sign(mean(sign(V(:,1))));
%tmp_stack = sign(V(1))*U;

tmp_stack = pws(indiv_tmp,true,false,1); tmp_stack = tmp_stack/norm(tmp_stack);
%tmp_stack = median(indiv_tmp,2,"omitnan"); tmp_stack = tmp_stack/norm(tmp_stack);
%tmp_stack = mean(indiv_tmp,2,"omitnan"); tmp_stack = tmp_stack/norm(tmp_stack);

winlen = size(indiv_tmp,1);
tdum = seconds((0:winlen-1)/Fs+shift);
%tdum = tdum - seconds(2^9);
li = length(indices);

for i = 1:li
    indiv_tmp(:,i) = ampFact*indiv_tmp(:,i)+i-1;
end

if stackFlag
    indiv_tmp = [indiv_tmp ampFact*tmp_stack+li];
end

figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
tiledlayout(1,1,"Padding","compact");
ax = nexttile();
plot(ax,tdum,indiv_tmp,'linewidth',2);
hold(ax,'on');
plot(ax,tdum,indiv_tmp(:,end),'k','linewidth',3);
ax.YDir = 'reverse';
axis tight
grid on;
xlabel('Seconds');
ylabel('Event Number');
tmp_stack = -(tmp_stack);
zoom on;