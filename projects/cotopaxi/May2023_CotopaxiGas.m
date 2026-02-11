%clear;
cd ~/igdata/SO2_COTOPAXI/;
DOAS_stations = ["2108126M1" "TAMD";...
    "D2J2230" "REFN";...
    "D2J2815" "REFS";...
    "D2J2835" "CAMI";...
    "I2J4969" "SANJ"];
%
lStations = length(DOAS_stations);
dirs = dir('2023.*');
lDirs = length(dirs);
TAMD = [];
REFN = TAMD;
REFS = TAMD;
CAMI = TAMD;
SANJ = TAMD;

%%
for i = 1:lDirs
    dir_ = dirs(i).name;
    cd(dir_);
    for j = 1:lStations
        code_ = DOAS_stations(j,1);
        kstnm_ = DOAS_stations(j,2);
        cd(code_);
        fName = strcat('PostFluxLog_',code_,'.txt'); try T_ = readtable(fName); catch; cd ..; continue; end;%,'FileType','text',"NumHeaderLines",2);
        if strcmp(kstnm_,"CAMI")
            CAMI = [CAMI; T_];
        elseif strcmp(kstnm_,"TAMD")
            TAMD = [TAMD; T_];
        elseif strcmp(kstnm_,"REFS")
            REFS = [REFS; T_];
        elseif strcmp(kstnm_,"REFN")
            REFN = [REFN; T_];
        elseif strcmp(kstnm_,"SANJ")
            SANJ = [SANJ; T_];
        end
        cd ../
    end
    cd ..
end

%%
close all;
tREFS = datetime(REFS.scandate)+REFS.scanstarttime;
tREFN = datetime(REFN.scandate)+REFN.scanstarttime;
tCAMI = datetime(CAMI.scandate)+CAMI.scanstarttime;
tSANJ = datetime(SANJ.scandate)+SANJ.scanstarttime;
tTAMD = datetime(TAMD.scandate)+TAMD.scanstarttime;

%ax__ = subplot(212);
figure();
semilogy(tREFS,REFS.flux__ton_day_,'o','LineWidth',2);
zoom on; grid on; hold on;
semilogy(tREFN,REFN.flux__ton_day_,'o','LineWidth',2);
semilogy(tCAMI,CAMI.flux__ton_day_,'o','LineWidth',2);
semilogy(tTAMD,TAMD.flux__ton_day_,'o','LineWidth',2);
semilogy(tSANJ,SANJ.flux__ton_day_,'o','LineWidth',2);
legend('refs','refn','cami','tamd','sanj','Location','northwest');
dumb_cbar = colorbar;
dumb_cbar.Visible = 'off';

