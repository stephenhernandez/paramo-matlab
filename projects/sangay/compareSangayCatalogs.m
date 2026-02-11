function ax = compareSangayCatalogs(dirname)
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Thursday 08 APR 2021

%%
if nargin < 1
    dirname = '~/research/now/sangay/';
end

%%
verdate = version('-date');
verdate = datenum(verdate);
versionFlag = false;
if verdate >= datenum(2017,09,14) %version 2017b
    % comment out any of these lines if they're giving you problems
    datetime.setDefaultFormats('default','dd-MMM-uuuu HH:mm:ss.SSS');
    format long g;
    
    set(groot,'defaultAxesFontSize',18);
    set(groot,'defaultColorbarFontSize',12);
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(groot,'defaultColorbarTickLabelInterpreter','latex');
    set(groot,'defaultTextarrowshapeInterpreter','latex');
    set(groot,'defaultTextboxshapeInterpreter','latex');
    set(groot,'defaultLineMarkerSize',12);
    set(groot,'defaultFigurePaperPositionMode','auto');
    set(groot,'defaultFigureRenderer','Painters');
    set(groot,'defaultFigureColor',[1 1 1]);
    versionFlag = true;
end

%%
if ~exist(dirname,'dir')
    fprintf(2,'directory %s does not exist, exiting',dirname);
    return;
end

%%
cd(dirname);

%% hugos infrasound cat
h_cat = load('hugo_sangay_catalog_clean.txt');
thugo = datenum(h_cat(:,1),h_cat(:,2),h_cat(:,3),h_cat(:,4),h_cat(:,5),h_cat(:,6));
if versionFlag
    thugo = datenum2datetime(thugo);
end
[thugo,hsi] = sort(thugo);
h_cat = h_cat(hsi,:);
h_cat = h_cat(:,7:11);

%% stephens seismic cat
s_cat = load('SangayAlternateCatalog.txt');
t = s_cat(:,1);
exmag = s_cat(:,3);
if versionFlag
    t = datenum2datetime(t);
end

%%
t_all = [t; thugo];
ampAll = [exmag; h_cat(:,1)];
i_all = [zeros(size(t)); ones(size(thugo))];

%%
[t_all,si] = sort(t_all);
i_all = i_all(si);
ampAll = ampAll(si);

%% HERE is the key part! Identify matches between both catalogs (change parameters if necessary)
mindiff = 115;
maxdiff = 125;
matchindex = 1;
matchI = seconds(diff(t_all)) >= mindiff & ...
    seconds(diff(t_all)) <= maxdiff & ...
    diff(i_all) == matchindex;
matchI2 = find(matchI) + 1;

%% find the actual index
refTime = min([min(thugo) min(t_all(matchI2))]);
lia = ismembertol(seconds(thugo - refTime),seconds(t_all(matchI2) - refTime),1,'DataScale',1);

%%
maxTI = max(thugo(lia));
hcatsize = size(h_cat,2);
ax = gobjects(hcatsize+1,1);
figh = gobjects(hcatsize,1);

%%
disp(sum(lia))
close all;
plotFlag = true;
if plotFlag
    ylabels = ["amplitude (infrasound)";"magnitud (seismic)";"azimuth";...
        "frequency";"velocity";"clt"];
    axn = 0;
    
    for i = 1:hcatsize
        figh(i) = figure('units','normalized','outerposition',[0 0 1 1]);
        axn = axn + 1;
        ax(axn) = axes;
        ax(axn).Parent = figh(i);
        if i ~= 1
            plot(ax(axn),thugo,h_cat(:,i),'.');
            hold on;
            plot(ax(axn),thugo(lia),h_cat(lia,i),'o');
            zoom on;
            ylabel(ylabels(axn));
            continue;
        end
        
        %%
        
        subplot(2,1,1,ax(axn));
        semilogy(ax(1),thugo,h_cat(:,i),'.');
        zoom on;
        hold on;
        semilogy(ax(1),thugo(lia),h_cat(lia,i),'o');
        ylabel(ylabels(axn));
        
        axn = axn + 1;
        ax(axn) = axes;
        ax(axn).Parent = figh(i);
        subplot(2,1,2,ax(axn));
        plot(t(t<=maxTI),exmag(t<=maxTI),'.');
        hold on;
        plot(t_all(matchI),ampAll(matchI),'o');
        ylabel(ylabels(axn));
        
    end
    linkaxes(ax,'x');
end
end

function tdt = datenum2datetime(tdn)
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
tdt = datetime(tdn,'ConvertFrom','datenum');
end
