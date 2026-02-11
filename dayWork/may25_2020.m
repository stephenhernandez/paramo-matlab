clear; close all; clc;
cd ~/research/now/pichincha_rrr/;
cd deformation/;
cd ggpa_raw/;

%%
load ggp_tilt.txt
t = datetime(ggp_tilt(:,1),ggp_tilt(:,2),ggp_tilt(:,3),ggp_tilt(:,4),ggp_tilt(:,5),ggp_tilt(:,6));

[t,tI] = sort(t);
rad = ggp_tilt(tI,7);
tang = ggp_tilt(tI,8);

%%
t = t(1:end-1);
Niir = 121;
Nmed = 121;
rad = cumsum(detrend(medfiltSH(zpkFilter(diff(rad),-inf,1/Niir,1,2,true),Nmed,true)));
tang = cumsum(detrend(medfiltSH(zpkFilter(diff(tang),-inf,1/Niir,1,2,true),Nmed,true)));

%%
n = 0;
years = (2015:2016)';
months = (01:12)';
for i = 1:length(years)
    year_ = years(i);
    for j = 1:length(months)
        n = n + 1;
        if n <= 18 && n >= 4
            month_ = months(j);
            
            %%
            fig = figure('units','normalized','outerposition',[0.15 0.15 0.75 0.75]);
            ax = subplot(211);
            plot(t,rad,'.-','linewidth',2);
            ylabel('Radial [$\mu r$]');
            grid on;
            title('Inclinometro GGPA');
            
            ax(2) = subplot(212);
            plot(t,tang,'.-','linewidth',2);
            ylabel('Tangencial [$\mu r$]');
            zoom on;
            grid on;
            
            %%
            linkaxes(ax,'x');
            if month_ == 12
                xlim([datetime(year_,month_,01) datetime(year_+1,01,01)]);
            else
                xlim([datetime(year_,month_,01) datetime(year_,month_+1,01)]);
            end
            cd ~/igreports/pichincha/tiltPlots/
            fname = strcat('GGPA_',datestr(datetime(year_,month_,01),'yyyy_mm'));
            print('-djpeg',fname);
            close all; 
        end
    end
end
