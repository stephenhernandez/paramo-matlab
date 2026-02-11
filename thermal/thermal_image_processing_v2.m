%for i in `ls *.fff`; do base=`echo ${i} |awk -F.fff '{print $1}'`; exiftool -b -RawThermalImage ${base}.fff > ${base}.tif; done
%
%https://github.com/VladimirSinitsin/flir_extractor/blob/main/fe_tools/fff_tools.py
%exiftool RV_rebeca_20180509_234410.fff -Emissivity -SubjectDistance ... <-- wildcard also optional, eg *.fff
% -AtmosphericTemperature -ReflectedApparentTemperature -IRWindowTemperature ...
% -IRWindowTransmission -RelativeHumidity -PlanckR1 -PlanckB -PlanckF ...
% -PlanckO -PlanckR2 -FieldOfView -j
%

clear; close all; clc;

%cd ~/igdata/fff_images/tif_images/
cd ~/data/reventador_thermal_images_different_vents/
load ReventadorBorderMain.mat
Dind = double(Dind);
Dind = Dind/norm(Dind);
[X,Y] = meshgrid((1:size(Dind,2))',(1:size(Dind,1))');
[Xq,Yq] = meshgrid((1:0.25:size(Dind,2))',(1:0.25:size(Dind,1))');
new = interp2(X,Y,Dind,Xq,Yq);

% close all;
% figure();
% imagesc(Xq(1,:)',Yq(:,1),new); zoom on; grid on; colorbar;

%
%cd ~/Desktop/IG-EPN_IRdata_Cotopaxivolcano_all/
%cd ~/Desktop/COTOPAXI_IR_Data_diario_raw_20240304/CTP_RUMHIR_20240304_0000-2357_fff

tic;
pr1 = 14615.877;    %PlanckR1
pr2 = 0.010745238;  %PlanckR2
pb = 1390.6; %PlanckB
po = -5913; %PlanckO
pf = 1; %PlanckF
%fov = 25; %not used
e = 0.97; %ice; 0.96; %Emissivity
r_temp = 0; %10; %20.0; %ReflectedApparentTemperature
raw_refl = pr1 / (pr2 * (exp(pb / (r_temp + 273.15)) - pf)) - po;

contrastFlag = true;
wienerFilter = true;
wienerKernel = [2,3];

files = dir('*.tif');
lFiles = length(files);
imData = NaN(240,320,lFiles);
imageMAD = NaN(lFiles,1);
t_image = NaT(lFiles,1);
columnPeak = NaN(lFiles,1);
columnTemp = columnPeak;
colMax = NaN(lFiles,1);
rowMax = colMax;

%
[X,Y] = meshgrid((1:320)',(1:240)');
[Xq,Yq] = meshgrid((1:0.25:320)',(1:0.25:240)');

%%
for i = 2019%1:lFiles
    fName = files(i).name;
    tStr = fName(11:25);
    t_image(i) = datetime(tStr,'InputFormat','yyyyMMdd_HHmmSS') + hours(5);
    imData = Tiff(fName,'r');
    imData = double(read(imData));
    imData = (imData - (1 - e) * raw_refl) / e;
    imData = pr1 ./ (pr2 .* (imData + po)) + pf;
    lnargI = imData > 0;

    imData = double(imData);
    imData(lnargI) = pb ./ log(imData(lnargI)) - 273.15;
    imData(~lnargI) = NaN;

    %%
    imageMAD(i) = mad(imData,1,'all');

    %%
    if wienerFilter
        imData = wiener2(imData,wienerKernel);
    end

    fullImage = imData;
    cutImage2 = fullImage(30:120,1:250);
    if contrastFlag
        windowMin = min(min(cutImage2));
        fullImage = fullImage - windowMin;
        cutImage2 = cutImage2 - windowMin;
        windowMax = max(max(cutImage2));
        cutImage2 = cutImage2./windowMax;
        fullImage = fullImage./windowMax;
    end

    fullImageResample = interp2(X,Y,fullImage,Xq,Yq);
    %BWs = edge(cutImage2,'sobel',0.01);
    BWs = edge(fullImageResample,'sobel',0.01);
    BWs = double(sign(BWs));
    BWs = BWs/norm(BWs);
    %cc = xcorr2(Dind,BWs);
    cc = xcorr2(new,BWs);

    [maxA,maxI] = max(cc);
    [~,colMax_] = max(maxA);
    rowMax_ = maxI(colMax_);
    rowMax(i) = rowMax_;
    colMax(i) = colMax_;

    %fullImageShifted = circshift(fullImage,[rowMax_-71 colMax_-250]);
    fullImageShifted = circshift(fullImageResample,[-(281-rowMax_) -(997-colMax_)]);
    fullImageShifted = fullImageShifted(500:end,280:end);
    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    imagesc(fullImageShifted);
    zoom on; grid on; colorbar;
end
toc;