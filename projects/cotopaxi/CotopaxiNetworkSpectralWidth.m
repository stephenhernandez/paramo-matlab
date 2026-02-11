clear; close all; %clc;

prewhiten = ~true;
newFs = 25.6;
shortWin = 20;  nOverlap1 = 4/5;
longWin = 640;  nOverlap2 = 29/32;
% newFs = 20;
% shortWin = 40;  nOverlap1 = 3/4;
% longWin = 600;  nOverlap2 = 9/10;
% newFs = 20;
% shortWin = 20;  nOverlap1 = 3/4;
% longWin = 600;  nOverlap2 = 9/10;

lfc = 1/(shortWin);
dW = 0; % window length (in Hz) for smoothing spectrum
detrendFlag = true;
gapFillValue = 0;
hfc = -inf; %floor(0.8*0.5*newFs);
npoles = 4;

% % kstnms = ["CO1V";"BREF";"BTAM";"BVC2";"BNAS";"COSE";"RUNA";"BMOR";"NAS2";"VC1";"TAMB"];
% kstnms = ["CO1V"; "BREF"; "BVC2"; "BTAM"; "BMOR"; "BNAS"; "NAS2"; "VC1"; ...
%     "COSE"; "SUCR"; "SLOR"; "TAMB"; "RUNA"; "BRRN"; "PITA"; "VCES"; "SRAM"];
% transferFlag = true;
% elementLats = [-0.665833;-0.66365;-0.66102;-0.67849;-0.72502;-0.6748;...
%     -0.6405;-0.642;-0.70834;-0.64005;-0.729811;-0.69366;-0.59431;-0.783167;...
%     -0.5585;-0.803284;-0.750222];
% elementLons = [-78.439167;-78.44082;-78.41417;-78.39913;-78.4591;-78.48751;...
%     -78.46767;-78.40312;-78.38206;-78.494483;-78.496745;-78.3575;-78.46706;...
%     -78.475389;-78.43267;-78.389817;-78.5669];
% elementElevs = [4930;4859;4408;4292;4337;3919;4022;4103;4162;3650;3600;4160;...
%     4000;3650;3710;4044;3073];
% dayVec = datetime(2023,05,12); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,["BHZ";"HHZ";"SHZ"],"EC","",true,true);

% kstnms = ["CTX1";"CTX2";"CTX3";"CTX4";"CTX5";"CTX6";"CTX7";"TAM1";"TAM2";"TAM3";"TAM4";"TAM5"]; transferFlag = false;
% elementLats = [-0.65787; -0.65682; -0.65602; -0.65567; -0.65543; -0.65421; -0.65452;-0.67766;-0.67942;-0.67846;-0.67788;-0.67970];
% elementLons = [-78.43968; -78.43797; -78.44049; -78.43839; -78.43647; -78.43969; -78.43798;-78.39572;-78.39528;-78.39729;-78.39800;-78.39720];
% elementElevs = [4643; 4625; 4581; 4587; 4551; 4542; 4534;4215; 4196; 4248; 4260; 4244];

kstnms = ["CTX1";"CTX2";"CTX3";"CTX4";"CTX5";"CTX6";"CTX7"]; transferFlag = false;
elementLats = [-0.65787; -0.65682; -0.65602; -0.65567; -0.65543; -0.65421; -0.65452];
elementLons = [-78.43968; -78.43797; -78.44049; -78.43839; -78.43647; -78.43969; -78.43798];
elementElevs = [4643; 4625; 4581; 4587; 4551; 4542; 4534];
% kstnms = ["TAM1";"TAM2";"TAM3";"TAM4";"TAM5"]; transferFlag = false;
% elementLats = [-0.67766;-0.67942;-0.67846;-0.67788;-0.67970];
% elementLons = [-78.39572;-78.39528;-78.39729;-78.39800;-78.39720];
% elementElevs = [4215; 4196; 4248; 4260; 4244];
dayVec = datetime(2023,05,12);
dayInc = 1;
C = loadWaveforms(dayVec,dayInc,kstnms,"DPZ","EC","",true,true,'~/data/nodes/CotopaxiNodeData2023');

% %kstnms = ["VCH1";"SN01";"SN02";"SN03";"SN04";"SN05";"SN06";"SN07";"SN08";"SN09";"SN10";"SN11";"SN12";"SN13";"SN14";"SN15";"SN16"]; transferFlag = ~true;
% kstnmsTrial = ["SN03";"SN04";"SN05";"SN07";"SN08";"SN11";"SN12";"SN13";"SN15"]; transferFlag = ~true;
% %dayVec = datetime(2018,06,26); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,"HHZ","9D","",true,true,'/Volumes/Samsung_T5/data/iguana/BROADBAND');
% dayVec = datetime(2018,07,26); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnmsTrial,["HHZ";"BHZ"],["9D";"EC"],"",true,true); 
% %dayVec = datetime(2023,05,11); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,["BHZ";"HHZ"],"EC","",true,true);
% kstnms = pull(C,'kstnm');
% [elementLats,elementLons,elementElevs] = metaDataFromStationList(kstnms);

% kstnms = "SAG1"; transferFlag = false;
% elementLats = [-2.193637; -2.193408; -2.193701; -2.193944; -2.193641];
% elementLons = [-78.099766; -78.099465; -78.099223; -78.099517; -78.099506];
% elementElevs = [1236; 1236; 1234; 1233; 1238];
% %dayVec = datetime(2022,01,15); dayInc = 2; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);
% dayVec = datetime(2024,07,09); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);
% %dayVec = datetime(2021,12,01); dayInc = 2; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);
% %dayVec = datetime(2021,12,25); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);
% %dayVec = datetime(2023,05,12); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);

% kstnms = "RUNA"; transferFlag = false;
% elementLats = [-0.5943264; -0.5944882; -0.5943; -0.5940137; -0.5942238];
% elementLons = [-78.4673559; -78.4670238; -78.4670167; -78.4670507; -78.4667395];
% elementElevs = [1236; 1236; 1234; 1233; 1238];
% dayVec = datetime(2023,05,11); dayInc = 2; C = loadWaveforms(dayVec,dayInc,kstnms,"BDF","EC",["01";"02";"03";"04";"05"],true,true);

% kstnms = ["BBIL";"BULB";"BPAT";"BMAS"]; transferFlag = true;
% elementLats = [-1.45004;-1.43713;-1.50731;-1.50281];
% elementLons = [-78.48759;-78.40965;-78.43727;-78.47942];
% elementElevs = [2696;3005;3768;2965];
% dayVec = datetime(2016,02,26); dayInc = 1; C = loadWaveforms(dayVec,dayInc,kstnms,"BHZ","EC","",true,true);

%%
%tDayStart = dateshift(C(1).ref,'start','day');
tBegs = pull(C,'ref');
tStops = tBegs + pull(C,'e');
tStart = dateshift(max(tBegs),'end','minute');
tEnd = dateshift(min(tStops),'start','minute');
%Ccut = cutWaveforms(C,tStart,0,tEnd-tStart);
Ccut = cutWaveforms(C,tStart+hours(12),0,hours(3));

%%
tw = shortWin*newFs*1;
Cf = detrendWaveforms(...
    syncWaveforms(...
    detrendWaveforms(...
    taperWaveforms(...
    interpolateWaveforms(...
    resampleWaveforms(...
    detrendWaveforms(...)
    Ccut),...
    newFs),gapFillValue),tw)),true,true));

%%
if transferFlag
    Cf = resampleWaveforms(Cf,newFs);
    Cf = padWaveforms(Cf);
    Cf = scaleWaveforms(transferWaveforms(nanGapWaveforms(Cf,gapFillValue),lfc,hfc,npoles,newFs,"vel",1,false),1e9);
else
    % Cf = detrendWaveforms(Cf);
    % Cf = filterWaveforms(Cf,-inf,floor(0.8*newFs/2));
    % Cf = downsampleWaveforms(Cf,10);
end
Cf = padWaveforms(Cf);

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
nfft = 2^(nextpow2(winlen)+3);
nfft2 = nfft/2;
w = blackmanharris(winlen); %
%w = tukeywin(winlen,2*(1-nOverlap1)); %
fshift = (-nfft/2:nfft/2-1)'*(newFs/nfft);
Ddiag = NaN(nfft,lT);
SW = Ddiag;
numCombos = lC + 0.5*lC*(lC-1);
CCpairs = NaN(nfft*1,lT,numCombos);
D3 = NaN(nfft,numCombos);
[Hbu,Fbu] = freqOperator(nfft,1,2,newFs,max([floor(npoles/2) 1])); %trial filter

%%
for i = 1:lT
    disp(i);
    d_ = squeeze(dCube(:,i,:));
    dCube2 = [];
    for kk = 1:lC
        d__ = d_(:,kk);
        if prewhiten
        d__ = tdNorm(d__,shortWin,1,newFs);
        d__ = fdWhiten(d__,lfc,hfc,dW,newFs,true,true);
        end
        dcut_ = cutWindows(d__,shortWin*newFs,nOverlap1,detrendFlag);
        dcut_ = dcut_.*w;
        D_ = fft(dcut_,nfft);
        dCube2 = cat(3,dCube2,D_);
    end

    %% now loop through each frequency
    lT2 = size(dCube2,2);
    for j = 1:nfft
        W = squeeze(dCube2(j,:,:));
        AverageNetworkCovarianceMatrix = W'*W;
        AverageNetworkCovarianceMatrix = AverageNetworkCovarianceMatrix/lT2;
        % [V,D] = eig(AverageNetworkCovarianceMatrix);
        % [Ddiag_,ind] = sort(diag(D),'descend');
        % V = V(:,ind);

        [V,Ddiag_] = svd(AverageNetworkCovarianceMatrix,"econ");
        Ddiag_ = diag(Ddiag_);
        SW_ = sum((0:lC-1)'.*Ddiag_)./sum(Ddiag_);
        SW(j,i) = SW_;
        Ddiag_ = Ddiag_(1);
        Ddiag(j,i) = Ddiag_;
        %CC = AverageNetworkCovarianceMatrix;
        V_ = V(:,1); CC = V_*V_'; %here we are filtering!!!!!!!!!!
        triuCC = triu(CC);
        triuCC = triuCC';
        triuCC = triuCC(:);
        D3(j,:) = triuCC(triuCC ~= 0)';
    end

    freqfilt = Hbu;
    CCpairs(:,i,:) = D3;
end
fOrig = fftshift(fshift,1);

%%
close all;
figure('units','normalized','outerposition',[0 0 1 1]);
clear ax__;
ax__ = subplot(211);
imagesc(ax__,(tMain),fshift,fftshift(Ddiag,1)); zoom on; grid on; colorbar; axis xy; ylim([0 newFs/2]);
set(ax__,'ColorScale','log');
ax__(2,1) = subplot(212);
imagesc(ax__(2),(tMain),fshift,(fftshift(SW,1))); zoom on; grid on; colorbar; axis xy; ylim([0 newFs/2]);
linkaxes(ax__,'xy');

lfc2 = 0.4;
hfc2 = lfc2*4;
npoles = 4;
sumSW = mean(SW(fOrig>=lfc2 & fOrig<=hfc2,:));
% figure(); plot(tMain,sumSW,'.'); zoom on; grid on;

tic;
%[DOA,AV,DT_DIFF] = network_covariance_array_processing(CCpairs,lfc2,hfc2,elementLats,elementLons,elementElevs,newFs,npoles);
[DOA,AV,DT] = network_covariance_array_processing(CCpairs,lfc2,hfc2,elementLats,elementLons,elementElevs,newFs,npoles);
%[DOA,AV,ELEVANGLE] = network_covariance_array_processing(CCpairs,lfc2,hfc2,elementLats,elementLons,elementElevs,newFs,npoles);
DOA(DOA>=360) = DOA(DOA>=360)-360;
%ELEVANGLE = 90 - ELEVANGLE;
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

figure('units','normalized','outerposition',[0 0 1 1]);
nStations = length(elementLats);
nCombos = 0.5*nStations*(nStations-1);
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
for i = 1:nStations
    kstnm1 = kstnms(i);
    for j = i+1:nStations
        kstnm2 = kstnms(j);
        n = n+1;
        ax___(n,1) = nexttile;
        dt = DT(n,:);
        plot(tMain,dt,'.');
        title(sprintf("%s-%s",kstnm1,kstnm2));
        hold on;
        plot([min(tMain) max(tMain)],[0 0],'k--');
    end
end
linkaxes(ax___,'x'); zoom on;

% figure();
% scatter(ELEVANGLE,DOA,[],sumSW,'filled'); zoom on; grid on;
%
% figure(); 
% SS = scatter(ELEVANGLE,AV,[],sumSW,'filled'); zoom on; grid on; set(gca,'YScale','log');
