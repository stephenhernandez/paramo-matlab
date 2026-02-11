% function [CCpairs,kstnms,elementLats,elementLons,elementElevs,newFs,shortWin,...
%     longWin,nOverlap1,nOverlap2,tMain,fshift,fOrig,SW,Ddiag] = ...
%     ArrayAndNetworkSVD(dayVec)
clear; close all; clc;
prewhiten = true;

% newFs = 20;
% shortWin = 30;  nOverlap1 = 1/2;
% longWin = 480;  nOverlap2 = 0;

% %Short Window seismic
% newFs = 25;
% shortWin = 6;  nOverlap1 = 5/6;
% longWin = 60;  nOverlap2 = 5/6;

% %Long Window seismic
% newFs = 25.6;
% shortWin = 20;  nOverlap1 = 4/5;
% longWin = 640;  nOverlap2 = 61/64;

%Seismic Array
newFs = 25.6;
shortWin = 20;
nOverlap1 = 4/5;
longWin = 640;
nOverlap2 = 29/32;

% %Infrasound Array 2
% newFs = 25.6;
% shortWin = 10;
% nOverlap1 = 4/5;
% longWin = 640;
% nOverlap2 = 29/32;

%Infrasound Array 1
% newFs = 40;
% shortWin = 1;
% nOverlap1 = 3/4;
% longWin = 64; %0;
% nOverlap2 = 29/32;

refEllipse = referenceEllipsoid('wgs84');
lfc = 1/(shortWin);
dW = 0; % window length (in Hz) for smoothing spectrum
detrendFlag = true;
gapFillValue = 0;
hfc = -inf; %floor(0.8*0.5*newFs);
npoles = 4;

% kstnms = ["CO1V";"BREF";"BTAM";"BVC2";"BNAS";"COSE";"RUNA";"BMOR";"NAS2";"VC1";"TAMB"];
% kstnms = ["CO1V"; "BREF"; "BVC2"; "BTAM"; "BMOR"; "BNAS"; "NAS2"; "VC1"; ...
%     "COSE"; "SUCR"; "SLOR"; "TAMB"; "RUNA"; "BRRN"; "PITA"; "VCES"; "SRAM";"TOMA"];
% transferFlag = false;
% elementLats = [-0.665833;-0.66365;-0.66102;-0.67849;-0.72502;-0.6748;...
%     -0.6405;-0.642;-0.70834;-0.64005;-0.729811;-0.69366;-0.59431;-0.783167;...
%     -0.5585;-0.803284;-0.750222;-0.498389];
% elementLons = [-78.439167;-78.44082;-78.41417;-78.39913;-78.4591;-78.48751;...
%     -78.46767;-78.40312;-78.38206;-78.494483;-78.496745;-78.3575;-78.46706;...
%     -78.475389;-78.43267;-78.389817;-78.5669;-78.434611];
% elementElevs = [4930;4859;4408;4292;4337;3919;4022;4103;4162;3650;3600;4160;...
%     4000;3650;3710;4044;3073;3300];
% kstnms = ["BREF"; "BVC2"; "BTAM"; "BMOR"; "BNAS"; "NAS2"; "VC1"; ...
%     "SUCR"; "SLOR"; "TAMB";"BRRN"; "PITA"; "VCES"; "SRAM";"TOMA"];
% transferFlag = false;
% elementLats = [-0.66365;-0.66102;-0.67849;-0.72502;-0.6748;...
%     -0.6405;-0.642;-0.64005;-0.729811;-0.69366;-0.783167;...
%     -0.5585;-0.803284;-0.750222;-0.498389];
% elementLons = [-78.44082;-78.41417;-78.39913;-78.4591;-78.48751;...
%     -78.46767;-78.40312;-78.494483;-78.496745;-78.3575;...
%     -78.475389;-78.43267;-78.389817;-78.5669;-78.434611];
% elementElevs = [4859;4408;4292;4337;3919;4022;4103;3650;3600;4160;...
%     3650;3710;4044;3073;3300];

% %good
% kstnms = ["CO1V"; "BREF"; "BVC2"; "BTAM"; "BMOR"; "BNAS"];
% transferFlag = false;
% dayVec = datetime(2023,05,11);
% dayInc = 2;
% C = loadWaveforms(dayVec,dayInc,kstnms,["BHZ";"HHZ";"SHZ"],"EC","",true,true);
% kstnms = pull(C,'kstnm');
% [elementLats,elementLons,elementElevs] = metaDataFromStationList(kstnms);

% kstnms = ["CHL1"; "CHL2"; "CERN"; "ECEN"; "ICHI"; "IPAN";"ICAN"];
% transferFlag = false;
% dayVec = datetime(2025,05,17);
% dayInc = 1;
% C = loadWaveforms(dayVec,dayInc,kstnms,["HHZ";"SHZ"],["EC";"OP"],["";"00"]);
% kstnms = pull(C,'kstnm');
% [elementLats,elementLons,elementElevs] = metaDataFromStationList(kstnms);

% % kstnms = ["CO1V";"BREF"];
% % transferFlag = false;
% % elementLats = [-0.665833;-0.66365];
% % elementLons = [-78.439167;-78.44082];
% % elementElevs = [4930;4859];
% % % kstnms = ["CO1V";"BVC2"];
% % % transferFlag = false;
% % % elementLats = [-0.665833;-0.66102];
% % % elementLons = [-78.439167;-78.41417];
% % % elementElevs = [4930;4408];
% % % kstnms = ["CO1V"; "BREF"; "BVC2"];
% % % transferFlag = false;
% % % elementLats = [-0.665833;-0.66365;-0.66102];
% % % elementLons = [-78.439167;-78.44082;-78.41417];
% % % elementElevs = [4930;4859;4408];
% % %dayVec = datetime(2023,05,12);

% kstnms = ["CO1V";"BREF";"BVC2";"BTAM";"BNAS";"NAS2";"VC1"; ...
%     "SUCR";"SLOR";"TAMB";"BRRN";"PITA";"VCES";"TOMA"];
% transferFlag = false;
% elementLats = [-0.665833;-0.66365;-0.66102;-0.67849;-0.6748;...
%     -0.6405;-0.642;-0.64005;-0.729811;-0.69366;-0.783167;...
%     -0.5585;-0.803284;-0.498389];
% elementLons = [-78.439167;-78.44082;-78.41417;-78.39913;-78.48751;...
%     -78.46767;-78.40312;-78.494483;-78.496745;-78.3575;...
%     -78.475389;-78.43267;-78.389817;-78.434611];
% elementElevs = [4930;4859;4408;4292;3919;4022;4103;3650;3600;4160;...
%     3650;3710;4044;3300];
% dayVec = datetime(2024,04,03);
% dayInc = 1;
% C = loadWaveforms(dayVec,dayInc,kstnms,["BHZ";"HHZ";"SHZ"],"EC","",true,true);

% kstnms = ["BBIL";"BULB";"BPAT";"BMAS";"BRUN"]; transferFlag = false;
% elementLats = [-1.45004;-1.43713;-1.50731;-1.50281;-1.41852];
% elementLons = [-78.48759;-78.40965;-78.43727;-78.47942;-78.42639];
% elementElevs = [2696;3005;3768;2965;2639];
% %dayVec = datetime(2016,02,26);
% %dayVec = datetime(2015,04,10);
% dayVec = datetime(2013,01,194);
% dayInc = 2;
% C = loadWaveforms(dayVec,dayInc,kstnms,"BHZ","EC","",true,true);
% %C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC","",true,true);

% kstnms = ["CTX1";"CTX2";"CTX3";"CTX4";"CTX5";"CTX6";"CTX7";"TAM1";"TAM2";"TAM3";"TAM4";"TAM5"]; transferFlag = false;
% elementLats = [-0.65787; -0.65682; -0.65602; -0.65567; -0.65543; -0.65421; -0.65452;-0.67766;-0.67942;-0.67846;-0.67788;-0.67970];
% elementLons = [-78.43968; -78.43797; -78.44049; -78.43839; -78.43647; -78.43969; -78.43798;-78.39572;-78.39528;-78.39729;-78.39800;-78.39720];
% elementElevs = [4643; 4625; 4581; 4587; 4551; 4542; 4534;4215; 4196; 4248; 4260; 4244];

% kstnms = ["VCH1";"SN02";"SN04";"SN05";"SN06";"SN11";"SN12";"SN13";"SN14"]; transferFlag = ~true;
% elementLats = [-0.7824;-0.804489;-0.85007;-0.85194;-0.83914;-0.802;-0.811195;-0.61969;-0.84368];
% elementLons = [-91.114;-91.08802;-91.12921;-91.14756;-91.16856;-91.12399;-91.13496;-91.07489;-91.15261];
% elementElevs = [1020;984;1073;1045;1008;935;957;10;913];
% % kstnms = ["TAM1";"TAM2";"TAM3";"TAM4";"TAM5"]; transferFlag = false;
% % elementLats = [-0.67766;-0.67942;-0.67846;-0.67788;-0.67970];
% % elementLons = [-78.39572;-78.39528;-78.39729;-78.39800;-78.39720];
% % elementElevs = [4215; 4196; 4248; 4260; 4244];
% dayVec = datetime(2018,09,01); dayInc = 1;
% C = loadWaveforms(dayVec,dayInc,kstnms,"HHZ",["EC";"9D"],"",true,true,'/Users/stephen/data/galapagos/iguana/BROADBAND/');
%
% kstnms = ["VCH1";"SN03";"SN04";"SN07";"SN08";"SN11";"SN12";"SN13";"SN14";"SN15";"SN16"]; transferFlag = ~true;
% elementLats = [-0.7824; -0.77429; -0.85007; -0.80719;-0.85279;-0.802;-0.811195;-0.61969;-0.84368;-0.88217;-0.69141];
% elementLons = [-91.114;-91.09002;-91.12921;-91.17005;-91.1926;-91.12399;-91.13496;-91.07489;-91.15261;-91.11234;-91.26646];
% elementElevs = [1020;707;1073;1067;661;935;957;10;913;500;17];
% % kstnms = ["TAM1";"TAM2";"TAM3";"TAM4";"TAM5"]; transferFlag = false;
% % elementLats = [-0.67766;-0.67942;-0.67846;-0.67788;-0.67970];
% % elementLons = [-78.39572;-78.39528;-78.39729;-78.39800;-78.39720];
% % elementElevs = [4215; 4196; 4248; 4260; 4244];
% dayVec = datetime(2018,09,01); dayInc = 1;
% C = loadWaveforms(dayVec,dayInc,kstnms,"HHZ",["EC";"9D"],"",true,true,'/Users/stephen/data/galapagos/iguana/BROADBAND/');

kstnms = ["CTX1";"CTX2";"CTX3";"CTX4";"CTX5";"CTX6";"CTX7"]; transferFlag = false;
elementLats = [-0.65787; -0.65682; -0.65602; -0.65567; -0.65543; -0.65421; -0.65452];
elementLons = [-78.43968; -78.43797; -78.44049; -78.43839; -78.43647; -78.43969; -78.43798];
elementElevs = [4643; 4625; 4581; 4587; 4551; 4542; 4534];
% kstnms = ["TAM1";"TAM2";"TAM3";"TAM4";"TAM5"]; transferFlag = false;
% elementLats = [-0.67766;-0.67942;-0.67846;-0.67788;-0.67970];
% elementLons = [-78.39572;-78.39528;-78.39729;-78.39800;-78.39720];
% elementElevs = [4215; 4196; 4248; 4260; 4244];
dayVec = datetime(2023,05,10);
dayInc = 1;
rawdatadir = "~/masa/backups/"; %'/Volumes/heliotropic/data/cotopaxi/CotopaxiNodeData2023';
C = loadWaveforms(dayVec,dayInc,kstnms,"DPZ","EC","",true,true,rawdatadir);

%kstnms = ["VCH1";"SN01";"SN02";"SN03";"SN04";"SN05";"SN06";"SN07";"SN08";"SN09";"SN10";"SN11";"SN12";"SN13";"SN14";"SN15";"SN16"]; transferFlag = false;
%dayVec = datetime(2018,06,26); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,"HHZ","9D","",true,true,'/Volumes/Samsung_T5/data/iguana/BROADBAND');
%dayVec = datetime(2023,05,11); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,["BHZ";"HHZ"],"EC","",true,true);

% kstnms = "SAG1"; transferFlag = false;
% elementLats = [-2.193637; -2.193408; -2.193701; -2.193944; -2.193641];
% elementLons = [-78.099766; -78.099465; -78.099223; -78.099517; -78.099506];
% elementElevs = [1236; 1236; 1234; 1233; 1238];
% %dayVec = datetime(2022,01,15); dayInc = 2; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);
% dayVec = datetime(2024,03,06); dayInc = 2; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);
% %dayVec = datetime(2021,12,01); dayInc = 2; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);
% %dayVec = datetime(2022,01,01); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);
% %dayVec = datetime(2021,12,25); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);
% %dayVec = datetime(2023,05,12); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);

% kstnms = "RUNA"; transferFlag = false;
% elementLats = [-0.5943264; -0.5944882; -0.5943; -0.5940137; -0.5942238];
% elementLons = [-78.4673559; -78.4670238; -78.4670167; -78.4670507; -78.4667395];
% elementElevs = [1236; 1236; 1234; 1233; 1238];
% %dayVec = datetime(2024,04,03);
% dayVec = datetime(2023,05,12);
% %dayVec = datetime(2023,05,11);
% %dayVec = datetime(2024,06,21);
% %dayVec = datetime(2024,03,05);
% %dayVec = datetime(2024,03,07);
% %dayVec = datetime(2024,02,29);
% dayInc = 1;
% C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);


tDayStart = dateshift(C(1).ref,'start','day');
tBegs = pull(C,'ref');
tStops = tBegs + pull(C,'e');
tStart = dateshift(max(tBegs),'end','minute'); %+seconds(20)-minutes(1)+hours(0);
%tEnd = tStart + hours(1); %
tEnd = dateshift(min(tStops),'start','minute');
Ccut = cutWaveforms(C,tStart,0,tEnd-tStart);
Ccut = syncWaveforms(Ccut,1,1,1);
%Ccut = cutWaveforms(C,tDayStart+hours(3),0,hours(1)); %tEnd-tStart);

%
tw = 2*newFs/lfc;
Cf = syncWaveforms(Ccut,false,true,true);
Cf = resampleWaveforms(Cf,newFs);
Cf = nanGapWaveforms(Cf,0);
Cf = padWaveforms(Cf);
% detrendWaveforms(...
%     (...
%     detrendWaveforms(...
%     taperWaveforms(...
%     interpolateWaveforms(...
%     resampleWaveforms(...
%     detrendWaveforms(...)
%     (Ccut)),...
%     newFs*1),gapFillValue),tw))));
% Cf = padWaveforms(Cf);

if ~prewhiten
    Cf = normalizeWaveforms(Cf,detrendFlag,false,true);
end

%%
if transferFlag
    Cf = resampleWaveforms(Cf,newFs);
    Cf = padWaveforms(Cf);
    Cf = nanGapWaveforms(Cf,gapFillValue);
    Cf = transferWaveforms(Cf,lfc,hfc,npoles,newFs,"vel",1,false);
    Cf = scaleWaveforms(Cf,1e9);
    Cf = padWaveforms(Cf);
else
    % Cf = detrendWaveforms(Cf);
    % Cf = filterWaveforms(Cf,-inf,floor(0.8*newFs/2));
    % Cf = downsampleWaveforms(Cf,10);
end


%% get data cube
lC = length(C);
dCube = [];
for i = 1:lC
    d_ = Cf(i).d;
    dcut_ = cutWindows(d_,longWin*newFs,nOverlap2,detrendFlag);
    dCube = cat(3,dCube,dcut_);
end

[~,~,endIndex] = cutWindows(d_,longWin*newFs,nOverlap2,detrendFlag);
tMain = getTimeVec(Cf(1));
tMain = tMain(endIndex);

winlen = shortWin*newFs;
[~,lT,lC] = size(dCube);
nfft = 2^(nextpow2(winlen)+1);
nfft2 = nfft/2;

%w = rectwin(winlen);
%w = kaiser(winlen,10);
w = blackmanharris(winlen); %
%w = tukeywin(winlen,2*(1-nOverlap1)); %

fshift = (-nfft/2:nfft/2-1)'*(newFs/nfft);
fOrig = fftshift(fshift,1);
Ddiag = NaN(nfft,lT);
SW = Ddiag;
nXpairs = 0.5*lC*(lC-1);
numCombos = lC + nXpairs;
CCpairs = NaN(nfft*1,lT,numCombos);
D3 = NaN(nfft,numCombos);

%%
for i = 1:lT %for each big window....
    disp(i);
    d_ = squeeze(dCube(:,i,:));
    dCube2 = [];
    for kk = 1:lC
        d__ = d_(:,kk);
        if prewhiten
            d__ = tdNorm(d__,shortWin,1,newFs);
            d__ = fdWhiten(d__,lfc,hfc,dW,newFs,true,true);
        end
        dcut_ = cutWindows(d__,shortWin*newFs,nOverlap1,~detrendFlag);
        dcut = normalizeWaveforms(dcut_,true);
        dcut_ = w.*dcut_;
        D_ = fft(dcut_,nfft);
        dCube2 = cat(3,dCube2,D_);
    end

    %% now loop through each frequency
    %dCube2 is the data from all frequecnies and all lC channels associated
    %with a SINGLE "long window"
    %dCube2: nRows = nfft, nPages = nChans, nColumns = nTimeSlices
    svdFlag = true;
    lT2 = size(dCube2,2);
    for j = 1:nfft
        W = squeeze(dCube2(j,:,:))'; %after squeeze/transpose: nRows = nChans, nColumns = no. short wins at FIXED freq.
        [Uorig,Ddiag_,~] = svd(W,"econ");
        Ddiag_ = diag(Ddiag_);
        %SW_ = 1-((2*Ddiag_(end)^2)/(Ddiag_(1)^2 + Ddiag_(2)^2));   %planarity
        SW_ = 1-(Ddiag_(3)^2/Ddiag_(1)^2);                          %linearity (v1)
        %SW_ = 1-(Ddiag_(end)^2/Ddiag_(1)^2);                       %linearity (v2)
        %SW_ = sum((0:lC-1)'.*Ddiag_)./sum(Ddiag_);
        SW(j,i) = SW_;
        Ddiag(j,i) = Ddiag_(1);
        U = Uorig(:,1); %<-- vector associated with largest singular value
        CC = U*U';
        triuCC = triu(CC)';
        diagtriuCC = diag(triuCC);
        triuCC = triuCC - diag(diagtriuCC);
        D3(j,1:nXpairs) = triuCC(triuCC ~= 0)';
        D3(j,nXpairs+1:end) = diagtriuCC';
    end
    CCpairs(:,i,:) = D3;
end

%%
close all;
figure('units','normalized','outerposition',[0 0 1 1]);
clear ax__;
ax__ = subplot(211);
imagesc(ax__,(tMain),fshift,fftshift(Ddiag,1));
zoom on; grid on; colorbar; axis xy; ylim([0 newFs/2]);
set(ax__,'ColorScale','log');
clim([0.5 2.5]);

ax__(2,1) = subplot(212);
imagesc(ax__(2),(tMain),fshift,fftshift(SW,1));
zoom on; grid on; colorbar; axis xy; ylim([0 newFs/2]);
clim([0.4 1]);
linkaxes(ax__,'xy');

% 4-12 for RUNA Infrasound
lfc2 = 0.4;
hfc2 = lfc2*4;
% lfc2 = 4;
% hfc2 = 12; %lfc2*4;
npoles = 4;
fI = fOrig>=lfc2 & fOrig<=hfc2;
sumSW = sum(SW(fI,:))/sum(fI);
% figure(); plot(tMain,sumSW,'.'); zoom on; grid on;

tic;
[DOA,AV,DT] = network_covariance_array_processing(CCpairs,lfc2,hfc2,elementLats,elementLons,elementElevs,newFs,npoles);
DOA(DOA>=360) = DOA(DOA>=360)-360;
toc;

clear ax_
figure('units','normalized','outerposition',[0 0 1 1]);
ax_ = subplot(211);
SS = scatter(tMain,DOA,[],sumSW,'filled');
zoom on; ax = gca; ax.YScale = 'lin'; zoom on;
grid on;
ax.Box = 'on';
c = colorbar;
c.Label.String = 'spectral width';
ylabel('direction of arrival');

ax_(2,1) = subplot(212);
SS = scatter(tMain,AV,[],sumSW,'filled'); zoom on; ax = gca; ax.YScale = 'log'; zoom on;
grid on;
ax.Box = 'on';
c = colorbar;
c.Label.String = 'spectral width';
ylabel('apparent velocity [m/s]');

% ax_(3,1) = subplot(313);
% SS = scatter(tMain,ELEVANGLE,[],sumSW,'filled'); zoom on; ax = gca; ax.YScale = 'lin'; zoom on;
% grid on;
% ax.Box = 'on';
% c = colorbar;
% c.Label.String = 'spectral width';
% ylabel('elevation angle');

linkaxes(ax_,'x');
sgtitle(sprintf("lfc: %f; hfc: %f; npoles: %d",lfc2,hfc2,npoles));

%%
figure('units','normalized','outerposition',[0 0 1 1]);
nStations = lC; %length(elementLats);
nCombos = nXpairs; %0.5*nStations*(nStations-1);
facts = factor(nCombos);
for i = 1:length(facts)
    ratio(i) = prod(facts(1:i))/prod(facts(i+1:end));
end

ratio = abs(ratio-1);
[~,minI] = min(ratio);
nrows = prod(facts(1:minI));
ncolumns = nCombos/nrows;
n = 0;
tiledlayout(nrows,ncolumns, 'Padding', 'compact', 'TileSpacing', 'compact');

kcmpnms = pull(Cf,'kcmpnm');
for i = 1:nStations
    if any(strcmp(kcmpnms,"BDF"))
        kstnm1 = strcat(kstnms(1),kcmpnms(i));
    else
        kstnm1 = kstnms(i);
    end

    for j = i+1:nStations
        if any(strcmp(kcmpnms,"BDF"))
            kstnm2 = strcat(kstnms(1),kcmpnms(j));
        else
            kstnm2 = kstnms(j);
        end
        n = n+1;
        ax___(n,1) = nexttile;
        dt = DT(n,:);
        % dt = medfiltSH(dt',10,false); %experimental
        % dt = flipud(dt);
        % dt = medfiltSH(dt,10,false); %experimental
        % dt = flipud(dt);
        plot(tMain,dt,'.');
        title(sprintf("%s-%s",kstnm1,kstnm2));
        hold on;
        plot([min(tMain) max(tMain)],[0 0],'k--');
    end
end
linkaxes(ax___,'x'); zoom on;

%% optional -- tries to get phase changes on a per frequency level to then
%% get time delays and slowness inversion (like Yardim said)
% d_ = distance(elementLats,elementLons,-0.683727,-78.436542,refEllipse);
% fIndex = 40;
% stationRef = 1;
% W = squeeze(dCube2(fIndex,:,:))';
% W = W./rssq(W);
% [V1,Ddiag_,~] = svd(W,'econ');
% V1 = V1(:,1);
% erg = V1'*W;
% %erg = erg./norm(erg);
% Wf = erg.*V1;
%
% AverageNetworkCovarianceMatrix = Wf*Wf';
% AverageNetworkCovarianceMatrix/lT2;
% [V2,D] = eig(AverageNetworkCovarianceMatrix);
% [~,ind] = sort(diag(D),'descend');
% V2 = V2(:,ind);
% V2 = V2(:,1);
% CC = V2*V2';    %here we are filtering!!!!!!!!!!
%
% timeDelay = unwrap(angle(CC(stationRef,:)'))./(2*pi*fOrig(fIndex));
% close all;
% figure();
% plot(d_-d_(stationRef),timeDelay,'.'); zoom on; grid on;
% figure();
% plot(d_-d_(stationRef),(d_-d_(stationRef))./timeDelay,'.'); zoom on; grid on;