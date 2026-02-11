function createTropomiFrames(galapagosFlag)
%create tropomi frames
dataHome = fullfile("~","masa","tropomi"); %tharp
cd(dataHome);

nowTime = dateshift(datetime('now')+hours(5),'start','day');
tStart = nowTime-8;
tEnd = nowTime-1;
days = (tStart:tEnd)';
ldays = length(days);

if galapagosFlag
    minLat = -4;
    maxLat = 3;
    minLon = -96;
    maxLon = -88;

    galapagosHome = fullfile(dataHome,'galapagos');
    if ~exist(galapagosHome,'dir')
        mkdir(galapagosHome);
    end

    if ~exist(fullfile(galapagosHome,'frames'),'dir')
        mkdir(fullfile(galapagosHome,'frames'));
    end
    load ~/igdata/denseGalapagosGrid.mat
else
    minLat = -6;
    maxLat = 3;
    minLon = -85;
    maxLon = -74;
    load ~/igdata/denseEcuadorGrid.mat
end

[Xq,Yq] = meshgrid((minLon+0.5:0.005:maxLon-0.5)',...
    (minLat+0.5:0.005:maxLat-0.5)');
%
SAcolor = 'k';
cmap = flipud(hot(512));
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');

%%
tDay = NaT(ldays,1);
for i = 1:ldays
    close all;
    day_ = days(i);
    tDay(i,1) = day_;
    fprintf("%s %d\n",day_,i);
    [yyyy,mm,dd] = datevec(day_);

    trialNames = [sprintf("S5P_DLR_NRTI_01_L3_SO2_%04d%02d%02d.nc",yyyy,mm,dd);...
        sprintf("S5P_DLR_NRTI_01_040201_L3_SO2_%04d%02d%02d.nc",yyyy,mm,dd);...
        sprintf("S5P_DLR_NRTI_01_040100_L3_SO2_%04d%02d%02d.nc",yyyy,mm,dd);...
        sprintf("S5P_DLR_OFFL_01_L3_SO2_%04d%02d%02d.nc",yyyy,mm,dd);...
        sprintf("S5P_DLR_OFFL_01_040201_L3_SO2_%04d%02d%02d.nc",yyyy,mm,dd);...
        sprintf("S5P_DLR_OFFL_01_040100_L3_SO2_%04d%02d%02d.nc",yyyy,mm,dd)];
    lTrialNames = length(trialNames);


    %%
    successFlag = false;
    nTry = 0;
    while ~successFlag & nTry < lTrialNames
        nTry = nTry + 1;
        fileName = trialNames(nTry);
        try
            dat = ncread(fileName,'sulfur_dioxide_total_column');
            successFlag = true;
            %break;
        catch
            try
                dat = ncread(fileName,'sulphur_dioxide_total_column');
                successFlag = true;
            catch
                continue;
            end
        end
        % % for j = 1:lTrialNames
        % %     fileName = trialNames(j);
        % %     fprintf("%s\n",fileName);
        % %     try
        % %         dat = ncread(fileName,'sulfur_dioxide_total_column');
        % %         break;
        % %     catch
        % %         try
        % %             dat = ncread(fileName,'sulphur_dioxide_total_column');
        % %         catch
        % %             continue;
        % %         end
        % %     end
        % % end

        %%
        dat = dat';

        lonG = ncread(fileName,'lon');
        latG = ncread(fileName,'lat');

        lonG = (sort(double(lonG)));
        latG = (sort(double(latG)));

        latI = latG >= minLat & latG <= maxLat;
        latG = (latG(latI));
        dat = dat(latI,:);

        lonI = lonG >= minLon & lonG <= maxLon;
        lonG = lonG(lonI);
        dat = dat(:,lonI);

        tCoverageEnd  = ncreadatt(fileName,'/','time_coverage_end');
        C = textscan(tCoverageEnd,"%04d-%02d-%02dT%02d:%02d:%02d.%f");
        dvec_ = cat(1,C{:});
        tCoverageEnd = datetime(dvec_(1:6)');

        [r,c] = meshgrid(lonG,latG);
        close all;
        figure('units','normalized','outerposition',[0 0 0.6 1]);

        minDat = min(min(dat,[],1,"omitnan"),[],2,"omitnan");
        dat = 1 + dat - minDat;
        dI = ~isfinite(dat);
        dat(dI) = 1;

        r = r'; c = c'; dat = dat';
        F = griddedInterpolant(r,c,dat,"spline","nearest");
        Xq = Xq';
        Yq = Yq';

        dat2 = F(Xq,Yq);
        dat2 = dat2 + minDat - 1;
        dI = dat2 <= 0 | ~isfinite(dat2);
        dat2(dI) = NaN;

        dat2 = medfilt2(dat2,[13 13]); %<-- optional
        SS = imagesc(Xq(:),Yq(:),dat2');
        axis xy; zoom on;
        cb = colorbar;
        clim([0 5]);

        Xq = Xq';
        Yq = Yq';

        axis equal;
        colormap(cmap);
        axis([minLon+0.5 maxLon-0.5 minLat+0.5 maxLat-0.5]);
        hold on;

        plot(lonEC,latEC,SAcolor,'linewidth',3);
        plot(lon_noec,lat_noec,'-','linewidth',3,'color',SAcolor);

        title(sprintf("%s",tCoverageEnd)); %2024

        xlabel('Lon.');
        ylabel('Lat.');
        set(gca,'ColorScale','log'); %2024
        clim([0.1 1]);

        % cotopaxi
        %-0.725719, -78.380584
        %-0.635769, -78.496492
        plot([-78.496492 -78.380584 -78.380584 -78.496492 -78.496492],...
            [-0.63 -0.63 -0.725719 -0.725719 -0.63],'m','linewidth',2);

        % sangay
        % -1.953487, -78.387435
        % -2.050068, -78.298177
        plot([-78.387435 -78.298177 -78.298177 -78.387435 -78.387435],...
            [-1.953487 -1.953487  -2.050068 -2.050068 -1.953487],'m','linewidth',2);

        % reventador
        % -0.036426, -77.689879
        % -0.122584, -77.603364
        plot([-77.689879 -77.603364 -77.603364 -77.689879 -77.689879],...
            [-0.036426 -0.036426  -0.122584 -0.122584 -0.036426],'m','linewidth',2);

        if galapagosFlag
            jpgName = fullfile(galapagosHome,'frames',char(strcat("frame_",datestr(tDay(i),"yyyy.mm.dd"),"_uncorrected.jpg")));
        else
            jpgName = fullfile(dataHome,'frames',char(strcat("frame_",datestr(tDay(i),"yyyy.mm.dd"),"_uncorrected.jpg")));
        end
        title(sprintf("%s",tCoverageEnd));
        print('-djpeg',jpgName);
        pause(2);

        if galapagosFlag
            saveName = fullfile(galapagosHome,'frames',char(strcat(fileName,'.upsampled.mat')));
            save(saveName,'Xq','Yq','dat2','tCoverageEnd');
        else
            saveName = fullfile(dataHome,char(strcat(fileName,'.upsampled.mat')));
            save(saveName,'Xq','Yq','dat2','tCoverageEnd');
        end
    end
end