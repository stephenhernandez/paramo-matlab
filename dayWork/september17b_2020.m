clear; close all; clc;
cd ~/research/now/sangay/
load('sangayWoodAndersonAmplitudes'); %data from regional sensors

kstnms = ["SAGA";...
    "BPAT";...
    "BMAS";...
    "BULB";...
    "BRUN";...
    "PUYO";...
    "TAMH";...
    "PORT";...
    "PKYU";...
    "TAIS"];

load SAGA_WoodAnderson; %SAGA data
WA = [z2pWA z2pWA_BPAT z2pWA_BMAS z2pWA_BULB z2pWA_BRUN z2pWA_PUYO z2pWA_TAMH...
    z2pWA_PORT z2pWA_PKYU z2pWA_TAIS];
RAW = [z2pRaw z2pWA_BPAT_raw z2pWA_BMAS_raw z2pWA_BULB_raw z2pWA_BRUN_raw...
    z2pWA_PUYO_raw z2pWA_TAMH_raw z2pWA_PORT_raw z2pWA_PKYU_raw z2pWA_TAIS_raw];
minAmp = 1e-3;
maxAmpSAGA = 2e1;
maxAmp = 2e-1;
minRaw = 0;
minRawSAGA = 7e3;
maxRawSAGA = 1.25e5;
waFlag = true;

%%
% clear; close all; clc;
% cd ~/research/now/sangay/;
% load('sangayDisplacementAmplitudes');
% kstnms = ["SAGA";...
%     "BPAT";...
%     "BMAS";...
%     "BULB";...
%     "BRUN";...
%     "PUYO";...
%     "TAMH";...
%     "PORT";...
%     "PKYU";...
%     "TAIS"];
% 
% %
% waFlag = true;
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_SAGA_raw),20/180),0.6,1.2,6);
% z2pWA_SAGA_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_SAGA_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_SAGA_raw),20/180));
% if ~waFlag
%     z2pWA_SAGA_raw = scaleWaveforms(z2pWA_SAGA_raw,1e9);
% end
% z2pWA_SAGA = max(abs(double(pull(z2pWA_SAGA_raw))))';
% z2pWA_SAGA_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_BPAT_raw),20/180),0.6,1.2,6);
% z2pWA_BPAT_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_BPAT_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_BPAT_raw),20/180));
% if ~waFlag
%     z2pWA_BPAT_raw = scaleWaveforms(z2pWA_BPAT_raw,1e9);
% end
% z2pWA_BPAT = max(abs(double(pull(z2pWA_BPAT_raw))))';
% z2pWA_BPAT_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_BMAS_raw),20/180),0.6,1.2,6);
% z2pWA_BMAS_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_BMAS_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_BMAS_raw),20/180));
% z2pWA_BMAS = max(abs(double(pull(z2pWA_BMAS_raw))))';
% z2pWA_BMAS_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_BULB_raw),20/180),0.6,1.2,6);
% z2pWA_BULB_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_BULB_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_BULB_raw),20/180));
% if ~waFlag
%     z2pWA_BULB_raw = scaleWaveforms(z2pWA_BULB_raw,1e9);
% end
% z2pWA_BULB = max(abs(double(pull(z2pWA_BULB_raw))))';
% z2pWA_BULB_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_BRUN_raw),20/180),0.6,1.2,6);
% z2pWA_BRUN_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_BRUN_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_BRUN_raw),20/180));
% if ~waFlag
%     z2pWA_BRUN_raw = scaleWaveforms(z2pWA_BRUN_raw,1e9);
% end
% z2pWA_BRUN = max(abs(double(pull(z2pWA_BRUN_raw))))';
% z2pWA_BRUN_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_PUYO_raw),20/180),0.6,1.2,6);
% z2pWA_PUYO_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_PUYO_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_PUYO_raw),20/180));
% if ~waFlag
%     z2pWA_PUYO_raw = scaleWaveforms(z2pWA_PUYO_raw,1e9);
% end
% z2pWA_PUYO = max(abs(double(pull(z2pWA_PUYO_raw))))';
% z2pWA_PUYO_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_TAMH_raw),20/180),0.6,1.2,6);
% z2pWA_TAMH_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_TAMH_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_TAMH_raw),20/180));
% if ~waFlag
%     z2pWA_TAMH_raw = scaleWaveforms(z2pWA_TAMH_raw,1e9);
% end
% z2pWA_TAMH = max(abs(double(pull(z2pWA_TAMH_raw))))';
% z2pWA_TAMH_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_PORT_raw),20/180),0.6,1.2,6);
% z2pWA_PORT_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_PORT_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_PORT_raw),20/180));
% if ~waFlag
%     z2pWA_PORT_raw = scaleWaveforms(z2pWA_PORT_raw,1e9);
% end
% z2pWA_PORT = max(abs(double(pull(z2pWA_PORT_raw))))';
% z2pWA_PORT_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_PKYU_raw),20/180),0.6,1.2,6);
% z2pWA_PKYU_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_PKYU_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_PKYU_raw),20/180));
% if ~waFlag
%     z2pWA_PKYU_raw = scaleWaveforms(z2pWA_PKYU_raw,1e9);
% end
% z2pWA_PKYU = max(abs(double(pull(z2pWA_PKYU_raw))))';
% z2pWA_PKYU_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% Sf = filterWaveforms(taperWaveforms(detrendWaveforms(z2pWA_TAIS_raw),20/180),0.6,1.2,6);
% z2pWA_TAIS_raw = transferWaveforms(intWaveforms(Sf),[],[],6,10,'disp',true,waFlag);
% z2pWA_TAIS_raw = differentiateWaveforms(taperWaveforms(detrendWaveforms(z2pWA_TAIS_raw),20/180));
% if ~waFlag
%     z2pWA_TAIS_raw = scaleWaveforms(z2pWA_TAIS_raw,1e9);
% end
% z2pWA_TAIS = max(abs(double(pull(z2pWA_TAIS_raw))))';
% z2pWA_TAIS_raw = max(abs(double(pull(Sf))))';
% clear Sf
% 
% WA = [z2pWA_SAGA z2pWA_BPAT z2pWA_BMAS z2pWA_BULB z2pWA_BRUN z2pWA_PUYO z2pWA_TAMH...
%     z2pWA_PORT z2pWA_PKYU z2pWA_TAIS];
% RAW = [z2pWA_SAGA_raw z2pWA_BPAT_raw z2pWA_BMAS_raw z2pWA_BULB_raw z2pWA_BRUN_raw...
%     z2pWA_PUYO_raw z2pWA_TAMH_raw z2pWA_PORT_raw z2pWA_PKYU_raw z2pWA_TAIS_raw];
% 
% if waFlag
%     minAmp = 1e-3;
%     maxAmpSAGA = 1e2;
%     maxAmp = 2e-1;
%     minRaw = 0;
%     minRawSAGA = 1e3;
%     maxRawSAGA = 2e5;
% else
%     minAmp = 1e0;       % nanometers
%     maxAmpSAGA = 2e5;   % nanometers
%     maxAmp = 1e3;       % nanometers
%     minRaw = 0;
%     minRawSAGA = 1e3;
%     maxRawSAGA = 2e5;
% end

%%
lK = size(WA,2);
b = NaN(2,lK-1);
sagab = NaN(2,2);

%
plotFlag = true;

%
close all; 
if plotFlag
    figure('units','normalized','outerposition',[0 0 0.8 1]);
end

for i = 1:lK
    wa = WA(:,i);
    raw = RAW(:,i);
    if i == 1
        gI1 = wa > minAmp & raw > minRawSAGA & raw < maxRawSAGA & wa <= maxAmpSAGA & isfinite(wa) & isfinite(raw);
        b_ = robustfit(log10(raw(gI1)),log10(wa(gI1)));
        sagab(1,1) = b_(1);
        sagab(2,1) = b_(2);
        sum(gI1)
        
        if ~plotFlag
            figure();
            loglog(raw(gI1),wa(gI1),'.'); zoom on; grid on; hold on;
            xsynth = logspace(floor(log10(min(raw(gI1)))),ceil(log10(max(raw(gI1)))),201)';
            ysynth = 10.^(b_(1) + b_(2).*log10(xsynth));
            hold on;
            pp = loglog(xsynth,ysynth,'-','linewidth',3,...
                'DisplayName',['$b_{1}: ',num2str(b_(1)),', b_{2}: ',num2str(b_(2)),'$']);
        end
        
        %         gI2 = wa > minAmp & raw > minRawSAGA & wa <= maxAmpSAGA & raw < maxRawSAGA & isfinite(wa) & isfinite(raw) & t6 >= datetime(2018,11,11);
        %         %b_ = robustfit(log10([raw(gI1); raw(gI2)]),log10([wa(gI1); wa(gI2)]));
        %         b_ = robustfit(log10(raw(gI2)),log10(wa(gI2)));
        %         sagab(1,2) = b_(1);
        %         sagab(2,2) = b_(2);
        %         loglog(raw(gI2),wa(gI2),'.'); zoom on; grid on;
        %
        %         if ~plotFlag
        %             xsynth = logspace(floor(log10(min(raw(gI2)))),ceil(log10(max(raw(gI2)))),201)';
        %             ysynth = 10.^(b_(1) + b_(2).*log10(xsynth));
        %             hold on;
        %             pp = loglog(xsynth,ysynth,'-','linewidth',3,...
        %                 'DisplayName',['$b_{1}: ',num2str(b_(1)),', b_{2}: ',num2str(b_(2)),'$']);
        %         end
        %
        %         sum(gI2)
        disp(sagab);
    else
        if strcmp(kstnms(i),"BULB")
            % during this time period, the digitizer gain is wrong in the
            % official response records of the IG, need to tell Andrea
            % Cordova.
            if ~exist('t6','var')
                t6 = refs;
            end
            bI = t6 >= datetime(2016,01,32) & t6 <= datetime(2017,01,91);
            raw(bI) = raw(bI)/4;
            wa(bI) = wa(bI)/4;
        end
        
        gI = wa > minAmp & raw > minRaw & wa <= maxAmp & isfinite(wa) & isfinite(raw);
        b_ = robustfit(log10(raw(gI)),log10(wa(gI)));
        
        %%
        b(1,i-1) = b_(1);
        b(2,i-1) = b_(2);
        if plotFlag
            subplot(3,3,i-1)
            loglog(raw(gI),wa(gI),'.','DisplayName',['data, N=',num2str(sum(gI))]); zoom on; grid on;
            xsynth = logspace(floor(log10(min(raw(gI)))),ceil(log10(max(raw(gI)))),201)';
            ysynth = 10.^(b_(1) + b_(2).*log10(xsynth));
            hold on;
            pp = loglog(xsynth,ysynth,'-','linewidth',3,...
                'DisplayName',['$b_{1}: ',num2str(b_(1)),', b_{2}: ',num2str(b_(2)),'$']);
            pp.Color(4) = 0.8;
            title(kstnms(i));
            legend('show','location','best');
            xlabel('raw [counts]');
            if waFlag
                ylabel('wood-anderson [mm]');
            else
                ylabel('displacement [nm.]');
            end
        end
    end
end

%%
if plotFlag
    suptitle('Robust (bisquare) fits: Sangay explosion amplitudes (counts) to equivalent Wood-Anderson amplitude (mm)');
end
