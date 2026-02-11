cd ~/research/now/cotopaxi/
load cotopaxi_template_times_for_extracting.mat
goodI = eqmag >= 1;

kcmpnms = ["BHZ";"BHN";"BHE";"HHZ";"HHN";"HHE";"SHZ"];
kstnms = ["BREF";"BNAS";"BTAM";"BVC2";"BMOR";"VCES";"VC1";"NAS2"];

lfc = 0.625;
hfc = 5;
finalFs = 50;
diffFlag = false;
A = [];
kstnm2= [];
kcmpnm2 = [];

for i = 1:length(kstnms)
    kstnm_ = kstnms(i);
    for j = 1:length(kcmpnms)
        kcmpnm_ = kcmpnms(j);
        try
            A_ = extractWaveforms(t(goodI)-seconds(1),seconds(25),kstnm_,kcmpnm_,...
                "EC","",true,true,[],1,true,[lfc,hfc,false,false,diffFlag]);
            A_ = resampleWaveforms(A_,finalFs);
            A = [A A_];
            kstnm2 = [kstnm2; kstnm_];
            kcmpnm2 = [kcmpnm2; kcmpnm_];
        catch
            continue;
        end
        %save('savedWaveforms_cotopaxi_template_times_for_extracting.mat');
    end
end
save("savedWaveforms_cotopaxi_template_times_for_extracting");

% clear; close all;
%
% tStart = datetime(2017,01,01);
% [eqType,evStatus,t,eqlat,eqlon,eqdepth,~,ids,stderr,azgap,nPhases,nMLv,...
%     timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
%     locMethod,earthModel,creationTime,agencyID,evMode,scbullmag] = readCat1();
% eqmag = scbullmag;
% horError = sqrt(eqlonerr.^2 + eqlaterr.^2);
%
% refEllipse = referenceEllipsoid('wgs84');
% d_ = distance(eqlat,eqlon,-0.683727,-78.436542,refEllipse)*1e-3;
% [stla,stlo,stel] = metaDataFromStationList("BREF");
% rI = t >= tStart & (contains(eqType,"lp",'ignorecase',true) | ...
%     contains(eqType,"vlp",'ignorecase',true) | ...
%     contains(eqType,"vt",'ignorecase',true) | ...
%     contains(eqType,"unk",'ignorecase',true) | ...
%     contains(eqType,"hb",'ignorecase',true)) ...
%     & d_ < 15 & eqmag >= 0.7 & eqmag < 3 & ...
%     stderr <= 1.5 & nPhases >= 4 & eqdepth <= 10 & nMLv >= 4 & ...
%     ~strcmp(evStatus,"preliminary") & azgap <= 270 & ...
%     (strcmp(earthModel,"cotopaxi") | strcmp(earthModel,"iasp91")) & ...
%     timerr <= 1.5 & eqmagerr <= 0.5 & horError <= 5 & eqdeptherr > 0;
% 
% figure(); plot(eqlon(rI),eqlat(rI),'.'); zoom on; grid on; axis equal;
% figure(); plot(t(rI),eqmag(rI),'.'); zoom on;
% 
% tGood = t(rI);
% aGood = 10.^eqmag(rI);
% tOrig = tGood;
% eqmag = log10(aGood);
% nDays = 7;
% 
% [rate,~,medianMagsFixedTimeWin] = ...
%     t2r(tOrig,days(nDays),eqmag);
% 
% figure();
% axx(1) = subplot(311);
% semilogy(tOrig,rate/nDays,'.'); grid on; hold on;
% axx(2) = subplot(312);
% semilogy(tOrig,medianMagsFixedTimeWin,'.'); grid on; hold on;
% axx(3) = subplot(313);
% semilogy(tOrig,(10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays,'.');
% zoom on; grid on; hold on;
% linkaxes(axx,'x');
% 
% %
% eqlat2 = eqlat(rI);
% eqlon2 = eqlon(rI);
% t2 = t(rI);
% id = ids(rI);
% eqmag2 = scbullmag(rI);
% lE = length(id);
% E = populateSeisCompStructure(lE);
% for i = 1:lE
%     disp(i);
%     id_ = id(i);
%     E_ = readSCBulletin(id_);
%     E(i) = E_;
% end
% 
% %
% newMags = pull(E,'mag');
% idGood = false(lE,1);
% brefTime = NaT(lE,1);
% brefAmp = NaN(lE,1);
% brefHypoDist = brefAmp;
% for i = 1:lE
%     E_ = E(i);
%     pPhases = E_.Pphases;
%     kstnms = pull(pPhases,"stnm");
%     times = pull(pPhases,"t");
%     [lia,locb] = ismember("BREF",kstnms);
%     if lia
%         MLv = E_.MLv;
%         brefTime(i) = times(locb);
%         kstnms = pull(MLv,"stnm");
%         amps = pull(MLv,"amp");
%         [lia,locb] = ismember("BREF",kstnms);
%         if lia
%             idGood(i) = true;
%             brefAmp(i) = amps(locb);
%             horDist = distance(stla,stlo,E_.lat,E_.lon,refEllipse)*1e-3;
%             brefHypoDist(i) = sqrt(E_.depth.^2 + horDist.^2);
%         end
%     end
% end
% 
% %
% goodI =~isnat(brefTime) & newMags > 0.5;
% figure();
% plot(brefTime(goodI),newMags(goodI),'.'); zoom on;
% 
% tGood = brefTime(goodI);
% aGood = 10.^newMags(goodI);
% tOrig = tGood;
% eqmag = log10(aGood);
% nDays = 7;
% 
% [rate,~,medianMagsFixedTimeWin] = ...
%     t2r(tOrig,days(nDays),eqmag);
% 
% figure();
% axx(1) = subplot(311);
% semilogy(tOrig,rate/nDays,'.'); grid on; hold on;
% axx(2) = subplot(312);
% semilogy(tOrig,medianMagsFixedTimeWin,'.'); grid on; hold on;
% axx(3) = subplot(313);
% semilogy(tOrig,(10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays,'.');
% zoom on; grid on; hold on;
% linkaxes(axx,'x');
% 
% %%
% diffFlag = false;
% finalFs = 50;
% A = extractWaveforms(tGood-seconds(5),seconds(20),"BREF",["BHZ";"BHN";"BHE"],...
%     "EC","",true,true,[],1,true,[2,8,false,false,diffFlag]);
% 
% Zf = A(:,1);
% Nf = A(:,2);
% Ef = A(:,3);
% 
% Zf = resampleWaveforms(Zf,finalFs);
% Nf = resampleWaveforms(Nf,finalFs);
% Ef = resampleWaveforms(Ef,finalFs);
% 
% %
% df = [pull(Zf); pull(Nf); pull(Ef)];
% Z2 = populateWaveforms(size(df,2));
% N2 = Z2;
% E2 = Z2;
% 
% for i = 1:size(df,2)
%     Z2(i) = dealHeader(Z2(i),df(1:1000,i),finalFs,tGood(i)-seconds(5));
%     N2(i) = dealHeader(N2(i),df(1001:2000,i),finalFs,tGood(i)-seconds(5));
%     E2(i) = dealHeader(E2(i),df(2001:3000,i),finalFs,tGood(i)-seconds(5));
% end
% T = [Z2 N2 E2]';
