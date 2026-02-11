function [fig,ax,S_,t] = helicorder(S,unit,lfc,hfc,rmsNorm,diffFlag,...
    detrendFlag,ampFactor,colorFlag,localTime,normFlag)
if nargin < 2; unit = 'hour'; end
if nargin < 3; lfc = -inf; end
if nargin < 4; hfc = -inf; end
if nargin < 5; rmsNorm = true; end
if nargin < 6; diffFlag = false; end
if nargin < 7; detrendFlag = true; end
if nargin < 8; ampFactor = []; end
if nargin < 9; colorFlag = []; end
if nargin < 10; localTime = false; end
if nargin < 11; normFlag = true; end

%%
if ~ismember(string(unit),["hour","minute","day"])
    fprintf(2,"inappropriate unit. must be day, hour, or minute\n");
    return;
end

%%
ampFlag = true;
if isempty(ampFactor)
    ampFlag = false;
end

if isempty(unit)
    unit = 'hour';
end

verboseFlag = false;

%%
S = differentiateWaveforms(S);
if any(isfinite([lfc hfc]))
    npoles = 4;
    S = detrendWaveforms(S);
    S = filterWaveforms(S,lfc,hfc,npoles);
    if ~isfinite(lfc) && detrendFlag
        if verboseFlag
            fprintf('detrending data\n');
        end
        S = detrendWaveforms(S);
    end
else
    if verboseFlag
        fprintf('no filtering requested\n');
    end
    if detrendFlag
        if verboseFlag
            fprintf('detrending data\n');
        end
        S = detrendWaveforms(S);
    end
end

if ~diffFlag
    S = intWaveforms(S);
end

%%
lS = size(S,1);
lw = 1;
for i = 1:lS
    S_ = S(i);

    %%
    if isnat(S_.ref) || ~S_.npts
        if verboseFlag
            fprintf('this trace empty, cant plot, continue\n');
        end
        continue;
    end

    t = getTimeVec(S_);
    if localTime
        t = t - hours(5);
        S_.ref = S_.ref - hours(5);
    end

    tStart = dateshift(t(1),'start',unit);
    tEnd = dateshift(t(end),'end',unit);

    S_ = cutWaveforms(S_,tStart,0,seconds(tEnd-tStart),false,true);
    Fs = round(1/S_.delta);
    t = getTimeVec(S_);

    if strcmp(unit,'hour')
        winlen = 3600*Fs;
        xStr = 'HH:MM:SS';
    elseif strcmp(unit,'minute')
        winlen = 60*Fs;
        xStr = 'segundos de minuto';
    else
        winlen = 86400*Fs;
        xStr = 'hora del dia';
    end

    %%
    nOverlap = 0;
    S_ = S_.d;
    nanI = isnan(S_);
    if rmsNorm
        normer = rms(S_(~nanI));
    else
        normer = max(S_(~nanI));
    end
    if ~normFlag
        normer = 1;
    end
    S_ = S_/normer;

    S_ = cutWindows(S_,winlen,nOverlap,false);
    totSlices = size(S_,2);
    if ~ampFlag
        totDuration = days(tEnd - tStart);
        ampFactor = totDuration*625./mean(sqrt(sum(S_.^2,"omitnan")),"omitnan")/totSlices;
        if diffFlag
            ampFactor = ampFactor/totDuration;
        end
    end
    
    %%
    nanI = logical(cutWindows(nanI,winlen,nOverlap,false));
    t = cutWindows(datenum(t),winlen,nOverlap,false);
    t = dn2dt(t);

    %%
    unitStarts = flipud(t(1,:)');
    totSlices = size(S_,2);
    fig = figure("units","normalized","outerposition",[0.1 0 0.8 1]);
    ax = gca;
    fig.Visible = 'off';
    hold(ax,'on');

    %%
    for j = 1:totSlices
        t_ = t(:,j);
        if diffFlag
            S__ = S_(:,j);
            if isempty(colorFlag)
                plot(t_(1:end-1) - t_(1),diff(S__)/Fs + totSlices - j + 1,'linewidth',lw);
            else
                plot(t_(1:end-1) - t_(1),diff(S__)/Fs + totSlices - j + 1,'linewidth',lw,'Color','w');
            end
            S_(:,j) = S__;
        else
            nans_ = nanI(:,j);
            S__ = S_(:,j);
            S__(nans_) = 0;
            S__ = ampFactor*detrend(S__)+totSlices - j + 1;
            S__(nans_) = NaN;
            if isempty(colorFlag)
                plot(t_ - t_(1),S__,'linewidth',lw);
            else
                plot(t_ - t_(1),S__,'linewidth',lw,'Color','w');
            end
            S_(:,j) = S__;
        end
    end

    %%
    ylim([0 totSlices+1]);
    zoom on;
    xlabel(ax,xStr);
    if strcmp(unit,'hour')
        xlim(ax,[seconds(0) hours(1)]);
    end
    titStr = {strcat(S.knetwk,'.',S.kstnm,'.',S.khole,'.',S.kcmpnm); ...
        ['zero-to-peak: ',num2str((S.depmax - S.depmin)/2)]};
    title(titStr);

    if totSlices > 1
        yticks = ax.YTick;
        yI = yticks > 0 & yticks <= length(unitStarts) & yticks == round(yticks);
        yticks = yticks(yI);
        ax.YTick = yticks;
        unitStarts(yticks);
        masterTimes = dateshift(unitStarts(yticks),'start','hour','nearest');
        tickLables = cellstr(datestr(masterTimes));

        zeroHours = hour(masterTimes) == 0;
        sumZeroHours = sum(zeroHours);
        if sumZeroHours
            findZeros = find(zeroHours);
            for ii_ = 1:length(sumZeroHours)
                boldLabel = tickLables(findZeros(ii_));
                boldLabel = ['\bf',char(boldLabel)];
                tickLables{findZeros(ii_)} = boldLabel;
            end
        end
        ax.YTickLabel = tickLables;
        ax(1).Box = 'on';
    end
end
