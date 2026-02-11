clear; close all; clc;
S = loadWaveforms(datetime(2021,12,01),2,"SAGA","HHZ");
Sf = detrendWaveforms(cutWaveforms(resampleWaveforms(intWaveforms(filterWaveforms(detrendWaveforms(differentiateWaveforms(S)),3/8,24)),100),dateshift(S(1).ref,'start','day')+hours(12),0,hours(24)));

d = Sf(1).d;
t = getTimeVec(Sf(1));
[dcut,startIndex,endIndex,badFlag,nwindows] = cutWindows(Sf.d,16*441*100/15,15/16,true);
close all; cd ~/animations/saga_drumbeats/;
!\rm -rf figure*

fig = figure('units','normalized','outerposition',[0 0 1 1]);
fig.Visible = 'off';
ax(1) = subplot(5,1,(1:3)'); ax(2) = subplot(5,1,(4:5)');

for i = 1:size(dcut,2)
    dcut_ = dcut(:,i); 

    si = startIndex(i); 
    ei = endIndex(i); 
    ei3 = round(15*length(dcut_)/16); 
    ei2 = si+ei3;
    h1 = plot(ax(1),t(1:ei2),d(1:ei2),'k'); 
    hold(ax(1),'on'); 
    h2 = plot(ax(1),t(ei2:end),d(ei2:end),'Color',ax(1).ColorOrder(1,:)); 
    axis(ax(1),'tight');
    
    ax(1).YTick = [];
    tcut = t(si:ei);
    h3 = plot(ax(2),tcut(1:ei3),dcut_(1:ei3),'k'); hold(ax(2),'on');
    h4 = plot(ax(2),tcut(ei3:end),dcut_(ei3:end),'Color',ax(2).ColorOrder(1,:)); axis(ax(2),'tight');
    p2p = 0.5*peak2peak(dcut_); h5 = sgtitle(sprintf('zero-to-peak: %g\n',p2p),'FontSize',18);
    ax(2).YTick = [];
    pause(1); fName = sprintf('figure_%04d.png',i);
    print('-dpng',fName);

    delete(h1); 
    delete(h2);
    delete(h3);
    delete(h4); 
    delete(h5);
end