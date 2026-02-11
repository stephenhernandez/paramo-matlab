% may22_2020
clear; close all; clc;
volcanoNames = "cotopaxi"; %["guagua","cotopaxi"];
catPlots = false;
plotFlag = false;
printFlag = false;
%for ii = 1:length(volcanoNames)
volcanoName = char(volcanoNames(1)); % 'guagua';
[tsheet,type,Tsp,amp,magnitude,codaDuration,period,refStation] = loadVolcanicEvents(volcanoName);

%% remove repeated matches
[tsheet,amp,removeIndices] = removeRepeatedMatches(tsheet,amp,30);
for kk = 1:length(removeIndices)
    rI = removeIndices{kk};
    type(rI) = [];
    Tsp(rI) = [];
    codaDuration(rI) = [];
    period(rI) = [];
    refStation(rI) = [];
end

%
if strcmp(volcanoName,'cotopaxi')
    tt = tsheet >= datetime(2022,01,01) & amp >= 1e2 & refStation == "BREF";
end

%
yearStart = (2017:2022)';
monthStarts = (1:12)';
mm = 0;

close all;
t2 = tsheet(tt);
type2 = type(tt);
uniqueTypes = unique(type2);
ltypes = length(uniqueTypes);
for k = 1:ltypes
    uniqueTypes_ = uniqueTypes(k);
    figure();
    nDays = 7;
    tI2 = strcmp(uniqueTypes_,type2);
    t3 = t2(tI2);
    plot(t2(tI2),t2r(t2(tI2),nDays)/nDays,'.');
    zoom on; grid on;
    sgtitle(string(uniqueTypes_));
    dayCount = days(dateshift(max(t3),'end','day') - dateshift(min(t3),'start','day'));
    fprintf("%s: ntot: %d, avg.: %f\n",uniqueTypes_,sum(tI2),sum(tI2)/dayCount);
end

%%
for i = 1:length(yearStart)
    year_ = yearStart(i);
    for j = 1:length(monthStarts)
        mm = mm + 1;
        VT(mm) = 0;
        LP(mm) = 0;
        month_ = monthStarts(j);
        tStart = datetime(year_,month_,01);
        tStartMonth(mm) = tStart;
        if month_ == 12
            tEnd = datetime(year_+1,01,01); % january of next year
            if catPlots
                if strcmp(volcanoName,'guagua')
                    regional_catalog_plots('pichincha',year_,month_,01,year_+1,01,01,'ALL',20,-30,70,'o','spanish',21,'default',50);
                else
                    regional_catalog_plots(char(volcanoName),year_,month_,01,year_+1,01,01,'ALL',20,-30,70,'o','spanish',21,'default',50);
                end
                close all;
            end
        else
            tEnd = datetime(year_,month_+1,01);
            if catPlots
                if strcmp(volcanoName,'guagua')
                    regional_catalog_plots('pichincha',year_,month_,01,year_,month_+1,01,'ALL',20,-30,70,'o','spanish',21,'default',50);
                else
                    regional_catalog_plots(char(volcanoName),year_,month_,01,year_,month_+1,01,'ALL',20,-30,70,'o','spanish',21,'default',50);
                end
                close all;
            end
        end

        tI = tsheet >= tStart & tsheet < tEnd & amp >= 8e2;
        disp(' ');
        disp('-----');
        %disp(strcat(datestr(tStart,'yyyy mmm'),': ',num2str(sum(tI))));
        [N_,edges] = histcounts(tsheet(tI),dateshift(tStart,'start','day'):dateshift(tEnd,'start','day'));
        disp(strcat(datestr(tStart,'yyyy mmm '),': ',num2str(sum(tI)),',mean:',num2str(nanmean(N_')),',ci:',num2str(prctile(N_',[5 95]))));

        %%
        uniqueTypes = unique(type(tI));
        ltypes = length(uniqueTypes);
        Ndays = days(tEnd - tStart);
        N = NaN(Ndays,ltypes);
        legString = [];
        for k = 1:ltypes
            uniqueTypes_ = uniqueTypes(k);
            legString = [legString; string(uniqueTypes_)];
            tI2 = tI & strcmp(uniqueTypes_,type);
            sumType = sum(tI2);
            [N_,edges] = histcounts(tsheet(tI2),dateshift(tStart,'start','day'):dateshift(tEnd,'start','day'));
            N(:,k) = N_';
            dayCount = days(dateshift(max(tsheet(tI2)),'end','day') - dateshift(min(tsheet(tI2)),'start','day'));
            disp(strcat(uniqueTypes_,': ',num2str(sumType),',mean:',num2str(nanmean(N_')),',ci:',num2str(prctile(N_',[5 95]))));
            if strcmp(uniqueTypes_,"VT")
                VT(mm) = mean(N_,"omitnan");
            elseif strcmp(uniqueTypes_,"LP")
                LP(mm) = mean(N_,"omitnan");
            end

        end
        disp('-----');
        disp(' ');

        %%
        if plotFlag
            cd ~/Desktop/;
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            ax = axes;
            bar(edges(1:end-1)',N,'stacked');
            zoom on;
            legend(legString,'Location','NorthWest');
            grid on;
            ylabel('Numero Diario');
            if printFlag
                imageParent = strcat('~/igreports/',volcanoName,'BarPlots/');
                existFlag = exist(imageParent,'dir');
                if ~existFlag
                    disp('Directory doesnt exist. Creating...');
                    mkdir(imageParent)
                end
                fname = strcat(imageParent,'yearPlot_',volcanoName,'_',datestr(tStart,'yyyy'),'_',datestr(tStart,'mm'));
                print('-djpeg',fname);
                close all;
            end
        end
    end
end

