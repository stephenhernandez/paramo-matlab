
if exist('s1','var')
    s2 = now;
    elapsedTime = 86400*(s2 - s1);

    for i = 1:2
        fprintf('\n');
    end
    fprintf('<strong>elapsed time: %f seconds</strong>\n',elapsedTime);
    Sold = S;
    snclOld = string(strcat(pull(Sold,'knetwk'),pull(Sold,'kstnm'),pull(Sold,'khole'),pull(Sold,'kcmpnm')));
end

%%
clearvars -except S Sold snclOld conn;
s1 = now;
[S,Supdated] = updateWaveforms(S,1200,true);
s2 = now;
elapsedTime = 86400*(s2 - s1);
nUpdated = length(Supdated);
fprintf('<strong>updated %d in %f seconds (%f events per second...)</strong>\n',nUpdated,elapsedTime,nUpdated/elapsedTime);

if nUpdated == 1
    if isnat(Supdated(1).ref)
        return;
    end
end

%%
tic;
newFs = 64;
lfc = 0.6;
hfc = 1.2;
%Sf = demeanWaveforms(S); toc; disp('demean');

tw = 0.01;
Sf = differentiateWaveforms(S); toc; disp('differentiate');
Sf = detrendWaveforms(Sf); toc; disp('detrend');
Sf = taperWaveforms(Sf,tw); toc; disp('taper');
Sf = filterWaveforms(Sf,lfc,hfc); toc; disp('filtered');
Sf = resampleWaveforms(Sf,newFs); toc; disp('resample');
Sf = detrendWaveforms(Sf); toc; disp('detrend');
Sf = intWaveforms(Sf); toc; disp('integrate');
Sf = detrendWaveforms(Sf); toc; disp('detrend');

%Sf = intWaveforms((resampleWaveforms(differentiateWaveforms(((S))),newFs))); toc;
snclMaster = string(strcat(pull(S,'knetwk'),pull(S,'kstnm'),pull(S,'khole'),pull(S,'kcmpnm')));
snclUpdate = string(strcat(pull(Supdated,'knetwk'),pull(Supdated,'kstnm'),pull(Supdated,'khole'),pull(Supdated,'kcmpnm')));
[lia,locb] = ismember(snclUpdate,snclMaster);

%%
kcmpnms = pull(S,'kcmpnm');
mySNCLs = [pull(S,'kstnm') kcmpnms pull(S,'knetwk') pull(S,'khole')]; toc;
updatedSNCLs = [pull(Supdated,'kstnm') pull(Supdated,'kcmpnm') pull(Supdated,'knetwk') pull(Supdated,'khole')]; toc;
npts = pull(Sf,'npts'); toc;
td = dateshift(dn2dt(now)+hours(5),'start','day');
toc;

%%
%noNotUpdated = true;
sumlia = sum(lia);
nsumlia = false;
legStr = [];
if sumlia
    I = false(size(npts));
    I(locb(lia)) = true;

    close all;
    figure(); hold on;
    semilogy(find(I),(seconds(npts(I)/newFs)),'o'); zoom on; grid on; toc;
    legStr = [legStr; string(['Updated: ',num2str(sumlia)])];

    if sumlia < length(snclMaster)
        %noNotUpdated = false;
        semilogy(find(~I),(seconds(npts(~I)/newFs)),'.'); zoom on; grid on; toc;
        legStr = [legStr; string(['Not Updated: ',num2str(sum(~I))])];
    end

    if exist('snclOld','var')
        lia_ = ismember(snclMaster,snclOld);
        nsumlia = sum(~lia_);
    end

    if nsumlia
        % new channels
        locb_ = find(~lia_);
        semilogy(locb_,(seconds(npts(locb_)/newFs)),'o'); zoom on; grid on; toc;
        %if
        legStr = [legStr; string(['New: ',num2str(nsumlia)])];
        for i = 1:length(locb_)
            fprintf(1,'Added: %s\n',snclMaster(locb_(i)));
        end
    end
    legend(legStr,'Location','NorthEastOutside');
end

%%
refTime = pull(Sf,'ref'); toc;
bTime = pull(Sf,'b'); toc;
durTime = pull(Sf,'e'); toc;
eTime = refTime + bTime + durTime; toc;
nowTime = dn2dt(now) + hours(5); toc;
latency = nowTime - eTime; toc;
lS = length(S); toc;

%%
% responseStructure = singleSNCLFreqResponse(mySNCLs,td,td+365,npts,newFs,'vel'); toc; disp('response info')
%
% %%
% Hd = zpkOperator(1,8,newFs,4);
% Hd = [Hd; zpkOperator(2,8,newFs,4)];
% Hd = [Hd; zpkOperator(4,8,newFs,4)];
%
% for j = 1:length(Hd)
%     Hd_ = Hd(j);
%     for i = 1:length(S)
%         dorig = Sf(i).d;
%         df = filter(Hd_,dorig);
%         Sf(i).d = df;
%     end
% end
% toc;

%%
updatedI = find(I);
sameI = find(~I);

plotFlag = true;
if plotFlag
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on; zoom on;
    for i = 1:lS
        ppOrig(i) = plot(eTime(i),i,'>');
    end
    for i = 1:lS
        if ismember(i,updatedI)
            plot(eTime(i),i,'>','MarkerFaceColor',ppOrig(i).Color);
            pp = plot(refTime(i),i,'k.');
        elseif ismember(i,sameI)
            pp = plot(refTime(i),i,'o','Color',[0.5 0.5 0.5]);
        elseif nsumlia
            if ismember(i,locb)
                pp = plot(refTime(i),i,'rs'); %,'Color',[0.5 0.5 0.5]);
            end
        end
        pp.Color(4) = 0.25;
    end

    for i = 1:lS
        pp = plot([refTime(i) eTime(i)],i*ones(2,1),'k-','linewidth',0.01);
        pp.Color(4) = 0.25;
    end
    grid on
    ax = gca;
    xLimits = ax.XLim;
    text(repmat(xLimits(1),lS,1),(1:lS)',strcat(pull(S,'kstnm'),repmat(",",lS,1),pull(S,'kcmpnm')));
    ax.YDir = 'reverse';
    toc;
end
