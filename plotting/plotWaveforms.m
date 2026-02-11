function [fig,ax] = plotWaveforms(varargin)
%
% plotWaveforms plot waveforms structure
%
% [ax,data,t] = plotWaveforms(varargin)
%
%[ax,data,t] = plotWaveforms(S,lfc,hfc,symbolStyle,npoles,secFlag,diffFlag,plotEnvFlag)
%[ax,data,t] = plotWaveforms(ax,S,lfc,hfc,symbolStyle,npoles,secFlag,diffFlag,plotEnvFlag)
%
% S = loadWaveforms(datetime(2020,01,14),01,"FER1","BHZ");
% S = filterWaveforms(resampleWaveforms(differentiateWaveforms(S),100),5,10);
% close all; locs = stalta(S,4,8,2,false,true);
%

%
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019
% Modified Tuesday, Jan 28, 2025

%%
nVarargin = length(varargin);
if nVarargin < 1
    disp('not enough input arguments');
    return;
end

%% parse inputs
firstArg = varargin{1};
axIn = isobject(firstArg);

functionDefaultsPrelim = {...
    populateWaveforms(),...     % S
    -inf,...                    % lfc
    -inf,...                    % hfc
    '-',...                     % symbolStyle
    4,...                       % npoles
    false,...                   % secFlag
    false,...                   % diffFlag
    false,...                   % plotEnvFlag
    'on'};                      % Visibility

%%
lDefaults = length(functionDefaultsPrelim);
if axIn
    functionDefaults = cell(1,lDefaults + 1);
    functionDefaults(2:end) = functionDefaultsPrelim;
    functionDefaults{1} = firstArg;
    optsToUse = functionDefaults;
    if nVarargin < 2
        disp('sorry, you need to supply data to plot');
        return;
    end
    optsToUse(1:nVarargin) = varargin;
    [ax,S,lfc,hfc,symbolStyle,npoles,secFlag,diffFlag,plotEnvFlag,Visibility] = deal(optsToUse{:});
else
    functionDefaults = functionDefaultsPrelim;
    optsToUse = functionDefaults;
    optsToUse(1:nVarargin) = varargin;
    [S,lfc,hfc,symbolStyle,npoles,secFlag,diffFlag,plotEnvFlag,Visibility] = deal(optsToUse{:});
end

%%
if isempty(lfc)
    lfc = -inf;
end

if isempty(hfc)
    hfc = -inf;
end

if isempty(symbolStyle)
    symbolStyle = '-';
end

if isempty(npoles)
    npoles = 4;
end

if isempty(secFlag)
    secFlag = false;
end

if isempty(diffFlag)
    diffFlag = false;
end

if isempty(plotEnvFlag)
    plotEnvFlag = false;
end

%%
if diffFlag
    S = differentiateWaveforms(S);
    S = detrendWaveforms(S);
end

%%
zeroPhaseFlag = false;
cornersfin = isfinite([lfc hfc]);
tStart = pull(S,'ref');
goodI = ~isnat(tStart);

%%
S = S(goodI);
lS = length(S);
if any(cornersfin)
    if cornersfin(1) && ~diffFlag
        S = detrendWaveforms(S);
        S = differentiateWaveforms(S);
    end

    %%
    S = filterWaveforms(S,lfc,hfc,npoles,[],zeroPhaseFlag);

    %%
    if cornersfin(1) && ~diffFlag
        S = intWaveforms(S);
    end
end

%%
if axIn
    if secFlag
        tStart = NaT;
        mints = seconds(0);
    else
        mints = NaT;
    end
    maxts = mints;

    legString = repmat("",lS,1);
    % plot over older axes
    for i = 1:lS
        hold(ax(i,1),'on')
        t = getTimeVec(S(i));
        disp(min(t));

        %%
        if S(i).gapFlag
            gapInfo = S(i).gapInfo;
            lg = size(gapInfo,1);
            if lg < 100
                S(i) = nanGapWaveforms(S(i),NaN,true);
            end
        end
        data = double(pull(S(i)));

        %%
        if secFlag
            t = t - t(1);
            mint = t(1);
            maxt = t(end);
            mints = min([mints mint]);
            maxts = max([maxts maxt]);
        else
            lims = ax(i).XLim;
            mint = t(1);
            maxt = t(end);
            mints = min([mints mint lims(1)]);
            maxts = max([maxts maxt lims(2)]);
        end

        if isempty(char(S(i).khole)) || strcmp(S(i).khole,"")
            legString(i) = strcat(S(i).knetwk,'.',S(i).kstnm,'.',S(i).kcmpnm);
        else
            legString(i) = strcat(S(i).knetwk,'.',S(i).kstnm,'.',S(i).khole,'.',S(i).kcmpnm);
        end

        if secFlag
            hp = plot(ax(i,1),seconds(t),data,symbolStyle,'linewidth',3,'DisplayName',legString(i));
        else
            gca = ax(i,1);
            %yyaxis(gca,'right');
            hp = plot(ax(i,1),t,data,symbolStyle,'linewidth',3,'DisplayName',legString(i));
        end
        %        legend(ax(i,1),oldString)
        fig = gcf;
        %legend(ax(i,1),'Visible','on','Box','off','AutoUpdate','off');

    end

    if ~secFlag
        lims = [mints maxts];
        for i = 1:lS
            ax(i,1).XLim = lims;
        end
    end
    return;
end
%%
if ~lS
    disp('no data to plot');
    return;
end

%%
fig = figure('units','normalized','outerposition',[0.08 0.08 0.84 0.84]);
fig.Visible = Visibility;
tiledlayout(lS,1,"Padding","compact","TileSpacing","none");
if secFlag
    keepXLim = [];
    if lS > 10
        nCols = ceil(lS/10);
        nRows = ceil(lS/nCols);
        tiledlayout(nRows,nCols, 'Padding', 'compact', 'TileSpacing', 'compact', 'TileIndexing', 'columnmajor');
        keepXLim = unique([keepXLim; nRows*(1:nCols)']);
    else
        keepXLim = lS;
    end
end

%%
ax = gobjects(lS,1);
legString = repmat("",lS,1);
if secFlag
    mints = seconds(0);
else
    mints = NaT;
end
maxts = mints;

n = 0;
for i = 1:lS
    n = n+1;
    ax(n,1) = nexttile;
    S_ = S(i);

    %%
    if S_.gapFlag
        gapInfo = S_.gapInfo;
        lg = size(gapInfo,1);
        if lg < 1000
            S_ = nanGapWaveforms(S_,NaN,true);
        end
    end

    %%
    t = getTimeVec(S_) - hours(0);

    %%
    data = pull(S_);

    %%
    tStart = t(1);
    tEnd = t(end);
    if secFlag
        t = t - tStart;
        mint = t(1);
        maxt = t(end);
        mints = min([mints mint]);
        maxts = max([maxts maxt]);
    else
        mints = min([mints tStart]);
        maxts = max([maxts tEnd]);
    end

    %%
    knetwk_ = S_.knetwk;
    kstnm_ = S_.kstnm;
    khole_ = S_.khole;
    kcmpnm_ = S_.kcmpnm;
    legString_ = sprintf("%s.%s.%s.%s",knetwk_,kstnm_,khole_,kcmpnm_);

    if secFlag
        legString_ = sprintf('%s, %s',string(legString_),sprintf('%s',tStart));
        hp = plot(ax(n,1),t,data,symbolStyle,'linewidth',1);
        hp.Color(4) = 0.9;
        ax(n).Title.String = legString_;
        ax(n).Title.FontSize = 12;
        xlim tight;

        if ~ismember(n,keepXLim)
            ax(n,1).XTickLabel = [];
        end
    else
        hp = plot(ax(n,1),t,data,symbolStyle,'linewidth',2,'DisplayName',legString_);
        legend(ax(n,1),'Visible','on','Box','off','AutoUpdate','off','FontSize',16);
        hp.Color(4) = 0.9;
        xlim tight;
        %legend(ax(n,1),legString_,'Visible','on','Box','off');
    end

    %%
    env = [];
    if plotEnvFlag
        hold(ax(n,1),'on');
        env = abs(hilbert(data));
        plot(ax(n,1),t,env,'k','linewidth',3);
    end
    ax(n,1).YRuler.Exponent = 0;
    legString(n) = legString_;
end

%%
linkaxes(ax,'x'); zoom on;
ax(end).XTickLabelRotation = 15;
ax(end).XLim = [mints maxts];

if ~secFlag
    if lS > 1
        for i = 1:lS
            if ~mod(i,2)
                ax(i).YAxisLocation = "right";
            end
            if i < lS
                ax(i,1).XTickLabel = [];
            end
        end
    end
end