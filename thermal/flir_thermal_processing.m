%for i in `ls *.fff`; do base=`echo ${i} |awk -F.fff '{print $1}'`; exiftool -b -RawThermalImage ${base}.fff > ${base}.tif; done
%
%https://github.com/VladimirSinitsin/flir_extractor/blob/main/fe_tools/fff_tools.py
%exiftool RV_rebeca_20180509_234410.fff -Emissivity -SubjectDistance ... <-- wildcard also optional, eg *.fff
% -AtmosphericTemperature -ReflectedApparentTemperature -IRWindowTemperature ...
% -IRWindowTransmission -RelativeHumidity -PlanckR1 -PlanckB -PlanckF ...
% -PlanckO -PlanckR2 -FieldOfView -j
%
%

clear;
close all;
tic;
%cd ~/igdata/fff_images/tif_images/
cd ~/Desktop/reventador_thermal_images_different_vents/
load ReventadorBorderMain.mat
Dind = double(Dind);
Dind = Dind/norm(Dind);

%cd ~/Desktop/IG-EPN_IRdata_Cotopaxivolcano_all/
%cd ~/Desktop/COTOPAXI_IR_Data_diario_raw_20240304/CTP_RUMHIR_20240304_0000-2357_fff
files = dir('*.tif');
lFiles = length(files);
imData = NaN(240,320,lFiles);
t_image = NaT(lFiles,1);
for i = 1:lFiles
    fName = files(i).name;
    tStr = fName(11:25);
    t_image(i) = datetime(tStr,'InputFormat','yyyyMMdd_HHmmSS');

    imageData = Tiff(fName,'r');
    imageData = read(imageData);
    imData(:,:,i) = double(imageData);
    %imData = cat(3,imData,imageData_);
end
t_image = t_image + hours(5);
toc;

%
pr1 = 14615.877;    %PlanckR1
pr2 = 0.010745238;  %PlanckR2
pb = 1390.6; %PlanckB
po = -5913; %PlanckO
pf = 1; %PlanckF
%fov = 25; %not used
e = 0.97; %ice; 0.96; %Emissivity
r_temp = 0; %10; %20.0; %ReflectedApparentTemperature

imDataCorr = imData;
raw_refl = pr1 / (pr2 * (exp(pb / (r_temp + 273.15)) - pf)) - po;
imDataCorr = (double(imDataCorr) - (1 - e) * raw_refl) / e;
imDataCorr = pr1 / (pr2 * (imDataCorr + po)) + pf;
lnargI = imDataCorr > 0;
sum(sum(sum(lnargI)))

imDataCorr = double(imDataCorr);
imDataCorr(lnargI) = pb ./ log(imDataCorr(lnargI)) - 273.15;
imDataCorr(~lnargI) = NaN;

imageMAD = NaN(size(imDataCorr,3),1);
for i = 1:size(imDataCorr,3)
    imageData = squeeze(imDataCorr(:,:,i));
    imageMAD(i) = mad(imageData,1,'all');
end

%
%%this is a test for reventador images, ignore in case of Coto
close all;
wienerFilter = true;
% imDataCorrFiltered = imDataCorr;
% for i = 1:11
%     fullImage = squeeze(imDataCorr(:,:,i));
%     if wienerFilter
%         fullImage = wiener2(fullImage,[2,3]);
%     end
%     imDataCorrFiltered(:,:,i) = fullImage;
% end
% refImage = median(imDataCorrFiltered,3);


columnPeak = NaN(lFiles,1);
columnTemp = columnPeak;
colMax = NaN(lFiles,1);
rowMax = colMax;
for i = 1:lFiles
    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    fullImage = squeeze(imDataCorr(:,:,i));
    if wienerFilter
        fullImage = wiener2(fullImage,[2,3]);
    end

    %fullImage = fullImage./refImage;
    t_image_ = t_image(i);


    cutImage2 = fullImage(30:120,1:250);
    cutImage2 = cutImage2 - min(min(cutImage2));
    cutImage2 = cutImage2./max(max(cutImage2));
    BWs = edge(cutImage2,'sobel',0.01);
    BWs = double(sign(BWs));
    BWs = BWs/norm(BWs);
    cc = xcorr2(Dind,BWs);

    [maxA,maxI] = max(cc);
    [~,colMax_] = max(maxA);
    rowMax_ = maxI(colMax_);
    rowMax(i) = rowMax_;
    colMax(i) = colMax_;

    fullImage = circshift(fullImage,[rowMax_-71 colMax_-250]);
    % if t_image_ <= datetime(2018,05,21)
    %     cutImage = fullImage(60:70,:);
    % elseif t_image_ >= datetime(2019,09,01)
    %     cutImage = fullImage(50:60,1:250);
    % else
    %     cutImage = fullImage(40:51,1:250);
    % end

    fullImage  = fullImage(31:110,21:220);
    cutImage = fullImage(20:60,:);
    ax(1) = subplot(8,1,1);
    meanTemp = mean(cutImage);
    plot(meanTemp,'.');
    [maxTemp,maxTempI] = max(meanTemp);
    columnPeak(i) = maxTempI;
    columnTemp(i) = maxTemp;

    hold on;
    plot(maxTempI,maxTemp,'p','linewidth',2);
    plot([maxTempI maxTempI],ax(1).YLim,'k--','LineWidth',2);
    %xlim([1 320])
    xlim([1 100])

    cdumb = colorbar;
    cdumb.Visible = 'off';
    grid on;
    titStr = sprintf('Column associated with Peak Temp: %d',maxTempI);
    title(titStr);
    ax(2) = subplot(8,1,(2:8)');
    imagesc(fullImage);
    zoom on; grid on; colorbar;

    clim([8 18]);
    colormap turbo;
    t_image_ = t_image(i);

    to = text(10,10,datestr(t_image_),'FontSize',20);
    to.BackgroundColor = 'w';
    linkaxes(ax,'x');
    imageName = sprintf('figure_%04d',i);
    print('-djpeg',imageName);
    pause(0.03);
    close all;
end
toc;

%%
close all
satFlag = ~true;
newZero = -38;
ranger = 8;
for i = 250:300%size(imDataCorr,3)
    imageData = squeeze(imDataCorr(:,:,i));

    if satFlag
        imageData = imageData - newZero;
        sI = imageData <= 0;
        imageData(sI) = 0;

        maxT = max(imageData(430:480,:),[],'all'); %imageData_(402);
        %maxT = max(max(imageData_));
        imageData = ranger*imageData/maxT;
    end

    %imageData_ = log10(imageData_);
    %imageData_ = edge(imageData_,'canny');
    %imageData_ = edge(imageData_,'sobel');
    %imageData_ = edge(imageData_,'Prewitt');
    %imageData_ = edge(imageData_,'Roberts',0.15);


    figure('units','normalized','outerposition',[0 0.2 0.5 1-0.2]);
    imagesc(imageData);
    zoom on; colorbar; axis equal; axis tight;

    colormap parula;
    %title(sprintf("max temp: %f",max(max(imageData_))));
    imDataCorr(:,:,i) = imageData;
    title(sprintf('%s\n',files(i).name));
end
figure(); imagesc(mean(imDataCorr,3,'omitnan')); colorbar; zoom on;

%%
% for i = 1:size(imDataCorr,3)
%     figure(i); axis tight;
%     fName = sprintf('frame_roberts015_%02d',i);
%     print('-djpeg',fName);
%     pause(2);
% end