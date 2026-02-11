%august30_2019
clear;
close all;
clc;

cd('/Users/stephen/research/now/dV/sierra_negra');
lfc = 0.2;
hfc = 0.8;
SNCL1 = {["SN02","HHZ","9D",""]}; %,["VCH1","HHZ","9D",""],["SN06","HHZ","9D",""],["SN07","HHZ","9D",""],["SN14","HHZ","9D",""],["SN04","HHZ","9D",""],["SN05","HHZ","9D",""],["SN08","HHZ","9D",""],["SN15","HHZ","9D",""],["SN16","HHZ","9D",""],["SN13","HHZ","9D",""]};
SNCL2 = {["SN03","HHZ","9D",""]}; %["SN06","HHZ","9D",""],["SN07","HHZ","9D",""],["SN14","HHZ","9D",""],["SN04","HHZ","9D",""],["SN05","HHZ","9D",""],["SN08","HHZ","9D",""],["SN15","HHZ","9D",""],["SN16","HHZ","9D",""],["SN13","HHZ","9D",""],["VCH1","HHZ","9D",""],["SN03","HHZ","9D",""]};
startTime = datetime(2018,04,01);
endTime = datetime(2018,08,31);
refTime = datetime(2018,06,26);

%%
for i = 1:length(SNCL1)
    sncl1_ = SNCL1{i};
    stnm1 = sncl1_(1);
    for j = 1:length(SNCL2)
        sncl2_ = SNCL2{j};
        stnm2 = sncl2_(1);
        
        %%
        dirName = char(strcat(stnm1,'_',stnm2));
        if ~exist(dirName,'dir')
            mkdir(dirName);
        end
        cd(dirName);
        
        %%
        [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = ...
            dV(startTime,endTime,refTime,sncl1_,sncl2_,lfc,hfc,0,20,0);
        y = medfiltSH(dVsymmetric,41,true);
        figure(); plot(t,y,'.'); zoom on;
        figure(); plot(dayStack(:,60)); hold on; plot(nanmean(dayStack,2),'linewidth',1);
        
        %%
        for k = 1:8
            fName = strcat('figure_',num2str(k));
            figure(k); 
            savefig(fName);
            print('-djpeg',fName);
        end
        save(dirName,'-v7.3')
        cd ../
        close all;
    end
end