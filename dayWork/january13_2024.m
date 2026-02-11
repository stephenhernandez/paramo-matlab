clear; close all; clc;
cd ~/data/SuperMouv2024/profil010/
%cd ~/data/SuperMouv2024/profil014/


PWD = pwd;
figureDir = 'figures';
if exist(fullfile(PWD,figureDir),'dir')
    disp('dir already exists');
else
    mkdir(figureDir);
end

%%
suFiles = dir('*.su');
%suFiles = suFiles(1350:1359); %for profile 14
lFiles = length(suFiles);
secDur = 5;
maxSeconds = 5;
FsOrig = 2000;
maxWinlen = 1 + FsOrig*secDur;
nTraces = 96;
bigMatrix = NaN(maxWinlen,nTraces,lFiles);


Nmed = 1;
tw = 40;
dx = 6.25; %spacing between hydrophones
dx2 = dx/2;
Dshot2hydrophone1 = dx; %this is an assumption, i have no clue if its right!

npoles = 4;
dumbTraces = 4;
cmap = [bone(256); sky(256)];
zeroPhaseFlag = false;
clipLevel = 2e0;
upsampleFactor = 1;


lfc1 = 8;
hfc1 = 256;

close all;
dtShot = 6; % 6 seconds between shots
traceToShift = 1; %variable, depends on ship velocity
nMidpoints = nTraces;
midpointIndex = (1:nMidpoints)';
lastMidpointIndex = max(midpointIndex);

for i = 1:lFiles%699%1:nTraces/2+1%lFiles
    disp(i);
    thisFile = char(suFiles(i).name);
    [Data,SegyTraceHeaders,SegyHeader]=ReadSu(thisFile,'endian','b');
    Data = Data(:,dumbTraces+1:dumbTraces+nTraces);

    distance_from_shot = Dshot2hydrophone1 + dx*(0:nTraces-1)';
    clear ax;
    Nmed2D = 1;

    % d2 = detrend(freqDiff(detrend(Data)));
    % d2 = taper(detrend(d2),tw);
    % d2 = zpkFilter(d2,-inf,1/10,1,4,false);

    %d2 = Data; %detrend(freqDiff(Data));
    %d2 = medfiltSH(d2,Nmed,zeroPhaseFlag);
    %d2 = taper(detrend(d2),tw);
    %d2 = normalizeWaveforms(d2);
    %d2 = (medfilt2(d2,[1 1]));
    %d2 = zpkFilter(d2,lfc1,hfc1,FsOrig,npoles,zeroPhaseFlag);
    %d2 = detrend(cumsum(detrend(d2)));
    %d2 = peak2rms(d2).*d2./fliplr(rms(d2));

    %
    close all;
    % d2 = resample(d2,upsampleFactor,1);
    % Fs = FsOrig*upsampleFactor;
    % tvec = (0:size(d2,1)-1)'/Fs;
    % fig1 = figure('units','normalized','outerposition',[0.2 0 0.4 1]);
    % fig1.Visible = 'on';
    % imagesc(distance_from_shot,tvec,d2); zoom on; grid on; c1 = colorbar;
    % colormap(cmap);
    % clim(clipLevel*[-1 1]);
    % ax(1) = gca; ax(1).XDir = 'reverse';
    % xlabel('distance [m]');
    % ylabel('twtt [sec.]');
    % title(sprintf('shot number: %d\n',SegyTraceHeaders(1).FieldRecord));

    %%
    %d2 = diff(Data); %detrend(freqDiff(detrend(Data)));
    d2 = Data; d2 = taper(detrend(d2),tw);
    %d2 = zpkFilter(d2,-inf,1/8,1,npoles,false);
    %d2 = zpkFilter(d2,-inf,1/4,1,npoles,false);
    d2 = zpkFilter(d2,10,-inf,FsOrig,npoles,0);

    % lfc2 = lfc1; hfc2 = hfc1;
    % d2 = freqDiff(Data,FsOrig);
    % d2 = medfiltSH(d2,Nmed,zeroPhaseFlag);
    % d2 = taper(detrend(d2),tw);
    % d2 = normalizeWaveforms(d2);
    % d2 = (medfilt2(d2,[Nmed2D Nmed2D]));
    % d2 = zpkFilter(d2,lfc2,hfc2,FsOrig,npoles,zeroPhaseFlag);
    % d2 = peak2rms(d2).*d2./fliplr(rms(d2));

    d2 = resample(d2,upsampleFactor,1);
    Fs = FsOrig*upsampleFactor;
    tvec = (0:size(d2,1)-1)'/Fs;

    experimentalFlag = false;
    if experimentalFlag
        d2 = normalizeWaveforms(d2);
        [data3,lags] = doAutoCorrFreqDom(d2); %doCrossCorrFreqDom(d2(:,1:end),repmat(d2(:,1),1,size(d2,2)),true);
        lI = lags >= 0 & lags <= 10000;
        data3 = normalizeWaveforms(data3(lI,:));

        [data4,lags] = doCrossCorrFreqDom(sign(data3(:,1:end)),repmat(data3(:,1),1,size(d2,2)),true);
        lI = lags >= 0 & lags <= 10000;
        d2 = normalizeWaveforms(data4(lI,:));
        d2 = peak2rms(d2).*d2./fliplr(rms(d2));

        fig2 = figure('units','normalized','outerposition',[0.6 0 0.4 1]);
        fig2.Visible = 'off';
        imagesc(distance_from_shot,tvec,d2); zoom on; grid on; c2 = colorbar;
        colormap(cmap);
        clim(10*clipLevel*[-1 1]);

        ax(2) = gca; ax(2).XDir = 'reverse';
        xlabel('distance [m.]');
        ylabel('twtt [sec.]');

        linkaxes(ax,'xy');
        ylim([0 maxSeconds]);
        title(sprintf('shot number: %d\n',SegyTraceHeaders(1).FieldRecord));

        % thisShot = thisFile(1:3);
        % fName = fullfile(figureDir,sprintf('figure_%s',thisShot));
        % print('-djpeg',fName);
    else
        % fig2 = figure('units','normalized','outerposition',[0.6 0 0.4 1]);
        % fig2.Visible = 'off';
        % imagesc(distance_from_shot,tvec,d2); zoom on; grid on; c2 = colorbar;
        % colormap(cmap);
        % clim(clipLevel*[-1 1]);
        %
        % ax(2) = gca; ax(2).XDir = 'reverse';
        % xlabel('distance [m.]');
        % ylabel('twtt [sec.]');
        %
        % %linkaxes(ax,'xy');
        % ylim([0 maxSeconds]);
        % title(sprintf('shot number: %d\n',SegyTraceHeaders(1).FieldRecord));
        %
        % % thisShot = thisFile(1:3);
        % % fName = fullfile(figureDir,sprintf('figure_%s',thisShot));
        % % print('-djpeg',fName);
    end

    winlen = size(d2,1);
    experiemtnFlag2 = true;
    if ~experiemtnFlag2
        bigMatrix(1:winlen,:,i) = d2(1:winlen,:);
        continue;
    end

    uPos = normalizeWaveforms(d2);
    sizeUPos = size(uPos);
    stf = normalizeWaveforms(uPos(70:200,1));
    tic;
    for jj = 1:nTraces
        u_ = normalizeWaveforms(uPos(:,jj));
        % [~,gpos_] = iterative_deconvolution(u_,stf,50,false); %flag2 true, slow
        % gpos3 = gpos_;
        % gI = gpos_ > 0;
        % gpos2 = gpos_(gI);
        % medamp_gpos = medfiltSH(gpos2,4,true);
        % gpos_(gI) = gpos_(gI)./medamp_gpos;
        % uPos_ = normalizeWaveforms(fftfilt(stf,gpos_));
        % uPos(:,jj) = uPos_;

        % uPos_ = iterative_deconvolution(u_,stf,50,false); %flag2 true, slow
        % uPos(:,jj) = uPos_;

        [~,gpos_] = iterative_deconvolution(normalizeWaveforms(uPos(:,jj)),stf,20,false); %flag2 false, fast
        gpos3 = gpos_;
        gI = gpos_ > 0;
        gpos2 = gpos_(gI);
        medamp_gpos = medfiltSH(gpos2,4,true);
        gpos_(gI) = gpos_(gI)./medamp_gpos;
        uPos_ = normalizeWaveforms(fftfilt(stf,gpos_));
        uPos(:,jj) = uPos_;

        %uPos_ = iterative_deconvolution(normalizeWaveforms(uPos(:,jj)),stf,100,false); %2min37sec
        %uPos(:,jj) = uPos_;
    end
    toc;
    % uPos_ = iterative_deconvolution(normalizeWaveforms(uPos(:)),stf,96*200,false);
    % uPos = reshape(uPos_,sizeUPos);
    bigMatrix(1:winlen,:,i) = uPos(1:winlen,:);
end

%%
bigOrig = bigMatrix;
clear bigMatrix;

%%
% nFold = 20;
% shifter = 5;
% shotProcessSequence = 1:lFiles-nFold+1;
% zot = zeros(winlen,shifter*length(shotProcessSequence));
% cumMidPoints = 0;
% for j = shotProcessSequence
%     disp(j); thisCluster = squeeze(bigOrig(1:winlen,:,j:j+nFold-1));
%     for jj = 1:nTraces
%         mask = ceil(190+(FsOrig*jj*dx/1490));
%         thisCluster(1:mask,jj,:) = 0;
%     end
%     for jj = 1:shifter
%         cumMidPoints = cumMidPoints + 1;
%         if jj == 1
%             CMP = zeros(winlen,nFold);
%             for i = 1:nFold
%                 CMP_ = circshift(thisCluster(:,:,i),[0 -((i-1)*shifter-(jj-1))]);
%                 CMP(:,i) = squeeze(CMP_(:,1));
%             end
%         else
%             CMP = zeros(winlen,nFold-1);
%             for i = 2:nFold
%                 CMP_ = circshift(thisCluster(:,:,i),[0 -((i-1)*shifter-(jj-1))]);
%                 CMP(:,i-1) = squeeze(CMP_(:,1));
%             end
%         end
%         cmp = normalizeWaveforms(CMP); ampFactor = 8e1;
%         cmp2 = resample(cmp,10,1);
%         raw_shifts = 0; previous_lag = 0;
%         for i = 2:size(cmp,2)
%             [cc,lags] = doCrossCorrFreqDom(cmp2(:,i-1),cmp2(:,i));
%             [~,maxI] = max(cc);
%             lags_ = lags(maxI)+previous_lag;
%             raw_shifts(i,1) = lags_;
%             previous_lag = lags_;
%         end
%         shifted_data = resample(apply_shifts(cmp2,raw_shifts),1,10);
%         zot(:,cumMidPoints) = pws(shifted_data,true,true,2); %tmp_stack = plot_family(shifted_data,1:size(shifted_data,2),30,10); %xlim(seconds([500 2000]));
%     end
% end
% figure(); imagesc(medfilt2(normalizeWaveforms(zot),[3 3]));
% zoom on; colorbar; clim(1e-5*clipLevel*[-1 1]); colormap bone;
