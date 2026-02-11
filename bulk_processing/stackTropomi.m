function stackTropomi()
dataHome = fullfile("~","masa","tropomi"); %tharp
cd(dataHome);
%
files = dir('*SO2*.mat');
lfiles = length(files);
minStack = zeros(1501,1601);
maxElem = 1501*1601;
saveFlag = true;

SE = strel("disk",7);
pivotFile = 269;

%
load("NoiseModelEcuadorTropomi");
outlieriness = outlieriness1;
outlierThreshold = 1.5;

%
LineWidth = 3;
shape_dir = fullfile("~","igdata","shape_files");
S = shaperead(fullfile(shape_dir,"volcanes_2.shp"));
SAcolor = 'k'; %[0.5 0.5 0.5];
cmap = flipud(hot(1024));
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');

smoothingWindow = 7;
figVisibility = 'off';
MeanDataForStacking = NaN(1501,1601,smoothingWindow);
conversionPolynomial = [1.30218179161098;...
    -3.62640591389595];

%
j = 1;
fileName = files(j).name;
load(fileName,'Xq','Yq');
t = NaT(lfiles,1);
tDay = t;

for j = 1:lfiles
    fileName = files(j).name;
    if j == 1
        load(fileName);
    else
        load(fileName,'dat2','tCoverageEnd');
    end

    tmpFname = char(fileName);
    if strcmp(string(tmpFname(1:4)),"S5P_")
        %fileName = strcat("S5P_DLR_NRTI_01_040100_L3_SO2_",yyyyStr,mmStr,ddStr,".nc");

        yyyy_ = str2double(tmpFname(end-10-14:end-7-14));
        month_ = str2double(tmpFname(end-6-14:end-5-14));
        day_ = str2double(tmpFname(end-4-14:end-3-14));
    else
        %fileName = char(strcat(yyyyStr,mmStr,ddStr,".S5P.TROPOMI.SO2.PGL.nc"));
        yyyy_ = str2double(fileName(1:4));
        month_ = str2double(fileName(5:6));
        day_ = str2double(fileName(7:8));
    end
    tDay(j) = datetime(yyyy_,month_,day_);
    t(j) = tCoverageEnd;
    disp(j);
end

%%
initialIndex = lfiles-2*smoothingWindow;
for j = initialIndex:lfiles
    fileName = files(j).name;
    if j == 1
        load(fileName);
    else
        load(fileName,'dat2','tCoverageEnd');
    end

    tmpFname = char(fileName);
    if strcmp(string(tmpFname(1:4)),"S5P_")
        %fileName = strcat("S5P_DLR_NRTI_01_040100_L3_SO2_",yyyyStr,mmStr,ddStr,".nc");
        yyyy_ = str2double(tmpFname(end-10-14:end-7-14));
        month_ = str2double(tmpFname(end-6-14:end-5-14));
        day_ = str2double(tmpFname(end-4-14:end-3-14));
    else
        %fileName = char(strcat(yyyyStr,mmStr,ddStr,".S5P.TROPOMI.SO2.PGL.nc"));
        yyyy_ = str2double(fileName(1:4));
        month_ = str2double(fileName(5:6));
        day_ = str2double(fileName(7:8));
    end
    tDay(j) = datetime(yyyy_,month_,day_);
    t(j) = tCoverageEnd;

    dat2 = dat2';
    badI = dat2 < 0 | ~isfinite(dat2);
    dat2(badI) = NaN;

    dat2 = medfilt2(dat2,[13 13]); %<-- optional! [11 11]
    badI2 = dat2 < 0 | ~isfinite(dat2);
    dat2(badI|badI2) = NaN;

    maxmax = sum(dat2,"all","omitnan");
    maxmax = maxmax./max([1 (sum(sum(~badI,1,"omitnan"),2,"omitnan")/maxElem)]);
    sumGas(j) = maxmax;

    dat_ = dat2;
    badI = dat_ <= 0 | ~isfinite(dat_);
    dat_(badI) = NaN;

    minVal = prctile(dat_(:),50);
    goodI = isfinite(dat_) & (dat_ > absFloor | dat_./outlieriness >= outlierThreshold);

    gasFilteredOut(j) = 10.^polyval(conversionPolynomial,log10(sum(dat_(~goodI),'all','omitnan')));
    dat_(~goodI) = NaN;
    maxmax = sum(dat_,"all","omitnan");
    gasFiltered(j) = maxmax;

    %
    fig = figure('units','normalized','outerposition',[0 0 0.6 1]);
    fig.Visible = figVisibility;
    SS = imagesc(Xq(:),Yq(:),dat_);
    axis xy; zoom on; colorbar;
    clim([0 10]);
    axis equal; colormap(cmap);
    axis([-83 -75 -5.5 2]);
    hold on;
    plot(lonEC,latEC,'k','linewidth',LineWidth);
    plot(lon_noec,lat_noec,'-','linewidth',LineWidth,'color',[0.5 0.5 0.5]);
    plot([S.X],[S.Y],'k--','linewidth',LineWidth-2);

    xlabel('Lon.');
    ylabel('Lat.');
    set(gca,'ColorScale','log');
    clim([0.1 1]);

    % cotopaxi
    %-0.725719, -78.380584
    %-0.635769, -78.496492
    plot([-78.496492 -78.380584 -78.380584 -78.496492 -78.496492],...
        [-0.63 -0.63 -0.725719 -0.725719 -0.63],'m','linewidth',LineWidth-1);

    % sangay
    % -1.953487, -78.387435
    % -2.050068, -78.298177
    plot([-78.387435 -78.298177 -78.298177 -78.387435 -78.387435],...
        [-1.953487 -1.953487  -2.050068 -2.050068 -1.953487],'m','linewidth',LineWidth-1);

    % reventador
    % -0.036426, -77.689879
    % -0.122584, -77.603364
    plot([-77.689879 -77.603364 -77.603364 -77.689879 -77.689879],...
        [-0.036426 -0.036426  -0.122584 -0.122584 -0.036426],'m','linewidth',LineWidth-1);

    gasCorrected_ = 10.^polyval(conversionPolynomial,log10(gasFiltered(j)));
    gasFiltered(j) = gasCorrected_;
    sgtitle(sprintf("%s, (%s); Total SO2: %s; Total Removed: %s",datestr(tDay(j)),datestr(t(j)),num2str(gasCorrected_,'%f'),num2str(gasFilteredOut(j))));

    PWD = pwd;
    jpgName = fullfile(PWD,'frames_filtered',strcat("frame_",datestr(tDay(j),"yyyy.mm.dd"),".jpg"));
    print('-djpeg',jpgName);
    pause(2);
    close all;

    if j <= initialIndex+smoothingWindow-1
        if smoothingWindow == 30
            MeanDataForStacking(:,:,j-initialIndex+1) = dat2;
        else
            MeanDataForStacking(:,:,j-initialIndex+1) = dat_;
        end
    else
        MeanDataForStacking_ = squeeze(MeanDataForStacking(:,:,2:end));
        MeanDataForStacking(:,:,1:smoothingWindow-1) = MeanDataForStacking_;
        if smoothingWindow == 30
            MeanDataForStacking(:,:,smoothingWindow) = dat2;
        else
            MeanDataForStacking(:,:,smoothingWindow) = dat_;
        end
    end

    if j < initialIndex+smoothingWindow-1
        gasCorrected1_ = 10.^polyval(conversionPolynomial,log10(sumGas(j)));
        fprintf("%s <strong>%f %f</strong>\n",t(j),gasCorrected1_,gasFiltered(j));
    else
        dat1_ = MeanDataForStacking;
        dat_ = sum(dat1_,3,"omitnan")/smoothingWindow; %mean daily value per pixel

        maxmax = sum(dat_,"all","omitnan");
        meanGasWeek(j) = maxmax;

        maxmax = sum(dat1_,"all","omitnan");
        meanGasWeek2(j) = maxmax/smoothingWindow;

        gasCorrected1_ = 10.^polyval(conversionPolynomial,log10(sumGas(j)));
        gasCorrected3_ = 10.^polyval(conversionPolynomial,log10(meanGasWeek(j)));
        fprintf("%s <strong>%f %f %f</strong>\n",t(j),gasCorrected1_,gasFiltered(j),gasCorrected3_);

        fig2 = figure('units','normalized','outerposition',[0 0 0.6 1]);
        fig2.Visible = figVisibility;
        SS = imagesc(Xq(:),Yq(:),dat_);
        axis xy; zoom on; colorbar;
        clim([0 10]);
        axis equal; colormap(cmap);
        axis([-83 -75 -5.5 2]);
        hold on;
        plot(lonEC,latEC,'k','linewidth',LineWidth);
        plot(lon_noec,lat_noec,'-','linewidth',LineWidth,'color',[0.5 0.5 0.5]);
        plot([S.X],[S.Y],'k--','linewidth',LineWidth-2);

        sgtitle(datestr(t(j)));
        xlabel('Lon.');
        ylabel('Lat.');
        set(gca,'ColorScale','log');
        clim([0.1 1]);

        % cotopaxi
        %-0.725719, -78.380584
        %-0.635769, -78.496492
        plot([-78.496492 -78.380584 -78.380584 -78.496492 -78.496492],...
            [-0.63 -0.63 -0.725719 -0.725719 -0.63],'m','linewidth',LineWidth-1);

        % sangay
        % -1.953487, -78.387435
        % -2.050068, -78.298177
        plot([-78.387435 -78.298177 -78.298177 -78.387435 -78.387435],...
            [-1.953487 -1.953487  -2.050068 -2.050068 -1.953487],'m','linewidth',LineWidth-1);

        % reventador
        % -0.036426, -77.689879
        % -0.122584, -77.603364
        plot([-77.689879 -77.603364 -77.603364 -77.689879 -77.689879],...
            [-0.036426 -0.036426  -0.122584 -0.122584 -0.036426],'m','linewidth',LineWidth-1);

        gasCorrected_ = 10.^polyval(conversionPolynomial,log10(meanGasWeek2(j)));
        sgtitle(sprintf("Mean over 7 days ending on: %s, Average Daily SO2: %s",datestr(tDay(j)),num2str(gasCorrected_,'%f')));

        PWD = pwd;
        if smoothingWindow == 7
            jpgName = fullfile(PWD,'frames_filtered_week',strcat("frame_",datestr(tDay(j),"yyyy.mm.dd"),"_mean7days.jpg"));
        elseif smoothingWindow == 30
            jpgName = fullfile(PWD,'frames_filtered_month',strcat("frame_",datestr(tDay(j),"yyyy.mm.dd"),"_mean7days.jpg"));
        end
        print('-djpeg',jpgName);
        pause(2);
        close all;
    end
end

%%
close all;
gasCorrected1 = 10.^polyval(conversionPolynomial,log10(sumGas));
gasCorrected3 = 10.^polyval(conversionPolynomial,log10(meanGasWeek));
gasCorrected4 = 10.^polyval(conversionPolynomial,log10(meanGasWeek2));
smoothWin = smoothingWindow;
medfiltTechnique = false;

%%
figure('units','normalized','outerposition',[0 0 1 1]);
clear ax;

semilogy(tDay,gasCorrected1,'o');
zoom on; grid on;
hold on;
semilogy(tDay,gasFiltered,'o');
hold on;
if medfiltTechnique
    plot(tDay,medfiltSH(gasCorrected1,smoothWin),...
        'linewidth',LineWidth,'Color',[0.5 0.5 0.5]);
    plot(tDay,medfiltSH(gasFiltered,smoothWin),'k','linewidth',LineWidth);
else
    plot(tDay,fftfilt(ones(smoothWin,1)/smoothWin,gasCorrected1),...
        'linewidth',LineWidth,'Color',[0.5 0.5 0.5]);
    plot(tDay,fftfilt(ones(smoothWin,1)/smoothWin,gasFiltered),'k','linewidth',LineWidth);
end

figure(1); hold on;
nDays = 30;
badGas = ~isfinite(gasCorrected3);
gasCorrected4(badGas) = 0;
gasSmoothed = fftfilt(ones(nDays,1)/nDays,gasCorrected4);
gasSmoothed(badGas) = NaN;
semilogy(tDay,gasSmoothed,'.'); zoom on; grid on;

if saveFlag
    clearvars -except outlieriness1 absFloor t tDay smoothingWindow gasCorrected* gasFiltered* meanGasWeek* sumGas*
    if smoothingWindow == 7
        save("NoiseModelEcuadorTropomi");
    elseif smoothingWindow == 30
        save("NoiseModelEcuadorTropomiMonthly");
    end
end