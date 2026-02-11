function [S,E] = writeStrongMotionRecordsNoPicks(eventID,noiseMinutes,...
    durMinutes,lfc,hfc,npoles,units)
feature('DefaultCharacterSet','ISO-8859-1'); %i dont know if this is necessary
if nargin < 1
    eventID = 'igepn2021xhpa.txt'; %28 nov 2021, M7.5 Peru, 10h50ish
end
if nargin < 2
    noiseMinutes = 1;
end

if nargin < 3
    durMinutes = 10;
end

if nargin < 4
    lfc = -inf;
end

if nargin < 5
    hfc = -inf;
end

if nargin < 6
    npoles = 4;
end

if nargin < 7
    units = 'WA';
end

%%
% clear; close all; clc;
% eventID = 'igepn2021xhpa.txt';
% %eventID = 'igepn2015qmro.txt';
% %eventID = 'igepn2016yvmv.txt';
% %eventID = 'igepn2013cuis.txt';
% eventID = 'igepn2019ghdt.txt';
% noiseMinutes = 1; durMinutes = 10; lfc = 0.5; hfc = 10; npoles = 4; units = 'WA';

%%
eventID = char(eventID);
eventID = eventID(1:13);
eventID = string(eventID);
E = readSCBulletin(string(eventID));
E = E(1);

%%
origdepth = E.depth;
origlat = E.lat;
origlon = E.lon;
origt = E.t;
origmag = E.mag;

%%
dayStart = dateshift(origt,'start','day');
yyyy = dayStart.Year;
mm = dayStart.Month;

%
yyyyStr = string(num2str(yyyy));
mmStr = sprintf("%02d",mm); %string(num2str(mm));

%%
finalFs = 100;
searchNets = ["EC";"IU";"OP";"II"];
fileList = existsToday(origt,true,true,[],searchNets);

%%
fileList = unique(fileList);
splitstring = split(fileList,".");
knetwks = splitstring(:,1);
kstnms = splitstring(:,2);
kholes = splitstring(:,3);
kcmpnms = splitstring(:,4);
kcmpnmBase = char(kcmpnms);
kcmpnmBase = string(kcmpnmBase(:,1:2));

%%
[stla,stlo,stel] = metaDataFromStationList(kstnms,knetwks,kcmpnms,kholes);
stel = stel*1e-3;
refEllipse = referenceEllipsoid('wgs84');
horDist = distance(origlat,origlon,stla,stlo,refEllipse)*1e-3;
d_ = sqrt(horDist.^2 + (origdepth+stel).^2);

minDist = 0; %10.^floor(log10(min(d_)));
maxDist = 2e3; %10.^ceil(log10(max(d_)));

% create filter
baseChanList = ["HN";"HH";"EN";"BH";"BL"];
goodI = d_ >= minDist & d_ <= maxDist & isfinite(stla) & ...
    ismember(kcmpnmBase,baseChanList) & ismember(knetwks,searchNets);

%%
% apply filter here
kstnms = kstnms(goodI);
kcmpnms = kcmpnms(goodI);
knetwks = knetwks(goodI);
kholes = kholes(goodI);
d_ = d_(goodI);
stla = stla(goodI);
stlo = stlo(goodI);
stel = stel(goodI);
horDist = horDist(goodI);

%%
[d_,goodI] = sort(d_);
kstnms = kstnms(goodI);
kcmpnms = kcmpnms(goodI);
knetwks = knetwks(goodI);
kholes = kholes(goodI);
stla = stla(goodI);
stlo = stlo(goodI);
stel = stel(goodI);
horDist = horDist(goodI);
kcmpnmBase = char(kcmpnms);
kcmpnmBase = string(kcmpnmBase(:,1:2)); %do not erase, i might need later to screen which files i write to sac...

%%
% recompute
[~,azs] = distance(origlat,origlon,stla,stlo,refEllipse);
[~,bazs] = distance(stla,stlo,origlat,origlon,refEllipse);

%
allSNCLs = strcat(knetwks,kstnms,kholes,kcmpnmBase);
[uniqSNCLs,ia] = unique(allSNCLs);
lUniqSNCLs = length(uniqSNCLs);
dScrambled = d_(ia);
[~,resortI] = sort(dScrambled);
uniqSNCLs = uniqSNCLs(resortI);

ncomps = 3; %length(componentsList);
S = populateWaveforms(ncomps*lUniqSNCLs);
n = 0;

%%
%saveRoot = fullfile('~','products','EcuadorStrongMotion'); %aki
saveRoot = fullfile("~","masa","strong_motion"); %tharp

origmag = round(origmag*100)/100;
eventDirName = sprintf("%s_M%3.2f_%s",eventID,origmag,datestr(origt,30));

writeDir = fullfile(saveRoot,yyyyStr,mmStr,eventDirName);
if ~exist(writeDir,'dir')
    mkdir(writeDir);
end

waveformsDir = fullfile(writeDir,strcat("waveforms_",upper(units)));
if ~exist(waveformsDir,'dir')
    mkdir(waveformsDir);
end

imageDir = fullfile(writeDir,strcat("images_",upper(units)));
if ~exist(imageDir,'dir')
    mkdir(imageDir);
end

%%
logFileName = fullfile(writeDir,strcat("log_",upper(units),".txt"));
logFileID = fopen(logFileName,'w');
summaryFileName = fullfile(writeDir,strcat("summary_",upper(units),".txt"));
summaryFileID = fopen(summaryFileName,'w');

%%
waFlag = false;
if strcmp(units,'acc')
    unitsStr = '$cm \cdot s^{-2}$';
    triggerThresh = 1/2; % cm/s^2
    titleFmt2 = "Distance: %3.2f km, Azimuth: %3.2f, Back Azimuth: %3.2f, Peak Component Acceleration: %3.2f (component: %s), Peak Vector Magnitude Acceleration: %3.2f, Est. Duration: %3.2f [sec.]";
    header1 = "Archivo de aceleracion, parte de la Red Nacional de Acelerografos (RENAC) Ecuador";
    fprintf(summaryFileID,'%s.%s.%s.%s %s %s %s %s %s %s %s %s \n',...
        'NET.STATION.LOCID.CMPNM MaxAcc MedianAcc TimeOfMaxAcc DIST AZ BAZ PeakVecACC EstimatedDuration');
    fprintf(summaryFileID,'\n');
elseif strcmp(units,'vel')
    unitsStr = '$cm \cdot s^{-1}$';
    triggerThresh = 0.1/2; % cm/s
    titleFmt2 = "Distance: %3.2f km, Azimuth: %3.2f, Back Azimuth: %3.2f, Peak Component Velocity: %3.2f (component: %s), Peak Vector Magnitude Velocity: %3.2f, Est. Duration: %3.2f [sec.]";
    header1 = "Archivo de velocidad, parte de la Red Nacional de Acelerografos (RENAC) Ecuador";
    fprintf(summaryFileID,'%s.%s.%s.%s %s %s %s %s %s %s %s %s \n',...
        'NET.STATION.LOCID.CMPNM MaxVel MedianVel TimeOfMaxVel DIST AZ BAZ PeakVecVel EstimatedDuration');
    fprintf(summaryFileID,'\n');
elseif strcmp(units,'disp')
    unitsStr = '$cm$';
    triggerThresh = 0.010/2; % cm
    titleFmt2 = "Distance: %3.2f km, Azimuth: %3.2f, Back Azimuth: %3.2f, Peak Component Displacement: %3.2f (component: %s), Peak Vector Magnitude Displacement: %3.2f, Est. Duration: %3.2f [sec.]";
    header1 = "Archivo de desplazamiento, parte de la Red Nacional de Acelerografos (RENAC) Ecuador";
    fprintf(summaryFileID,'%s.%s.%s.%s %s %s %s %s %s %s %s %s \n',...
        'NET.STATION.LOCID.CMPNM MaxDisp MedianDisp TimeOfMaxDisp DIST AZ BAZ PeakVecDisp EstimatedDuration');
    fprintf(summaryFileID,'\n');
elseif strcmp(units,'WA')
    waFlag = true;
    unitsStr = '$mm$';
    triggerThresh = 0.03; % 0.03 mm is a thresh for WA...?%0.010/2; % cm
    titleFmt2 = "Distance: %3.2f km, Azimuth: %3.2f, Back Azimuth: %3.2f, Peak Component Displacement: %3.2f (component: %s), Peak Vector Magnitude Displacement: %3.2f, Est. Duration: %3.2f [sec.]";
    header1 = "Archivo de desplazamiento, parte de la Red Nacional de Acelerografos (RENAC) Ecuador";
    fprintf(summaryFileID,'%s.%s.%s.%s %s %s %s %s %s %s %s %s \n',...
        'NET.STATION.LOCID.CMPNM MaxWADisp MedianWADisp TimeOfMaxDisp DIST AZ BAZ PeakVecDisp EstimatedDuration');
    fprintf(summaryFileID,'\n');
end

%%
headerVersion = 1;
fontSize = 14;
minDur = 1;
Visibility = 'off';

%%
peakAmps = NaN(lUniqSNCLs,1);
estDurs = peakAmps;
individualSNCLsProcessed = 0;

%%
for i = 1:lUniqSNCLs
    ncomps = 3; %reset, do not erase
    k = 0;
    thisSNCL = uniqSNCLs(i);
    lia = ismember(allSNCLs,thisSNCL);
    if sum(lia) ~= ncomps
        fprintf(1,"something went wrong with: %s, investigate!\n",thisSNCL);
        continue;
    end
    locb = find(lia);
    S3 = populateWaveforms(ncomps);
    for j = 1:ncomps
        tic;
        thisi = locb(j);
        knetwks_ = knetwks(thisi);
        kstnms_ = kstnms(thisi);
        kholes_ = kholes(thisi);
        kcmpnms_ = kcmpnms(thisi);
        sncldist = d_(thisi);
        snclaz = azs(thisi);
        snclbaz = bazs(thisi);
        nextDay = dateshift(origt,'end','day');
        dayInc = 1;
        MAXDUR = 20;
        if nextDay-origt < minutes(MAXDUR)
            %orig time is too close to next day, read two days
            dayInc = 2;
        end
        S1 = loadWaveforms(dayStart,dayInc,kstnms_,kcmpnms_,...
            knetwks_,kholes_,true,true);

        if isnat(S1.ref)
            fprintf('skipping %s\n',kstnms_);
            elapsed_time = toc;
            fprintf(logFileID,'Unable to load: %s.%s.%s.%s, exec. time: %d [seconds]\n',...
                knetwks_,kstnms_,kholes_,kcmpnms_,elapsed_time);
            continue;
        end

        n = n + 1;
        S1 = detrendWaveforms(S1);
        S1 = nanGapWaveforms(S1,0);
        S1 = cutWaveforms(S1,origt-minutes(noiseMinutes),0,minutes(MAXDUR));
        S1.evla = origlat;
        S1.evlo = origlon;
        S1.evdp = origdepth;
        S1.dist = sncldist;
        S1.az = snclaz;
        S1.baz = snclbaz;
        S1.gcarc = km2deg(sncldist);

        %
        S1 = taperWaveforms(S1,ceil(2*finalFs/lfc));
        if sum(~isfinite(S1.d))
            n = n-1;
            elapsed_time = toc;
            fprintf(logFileID,'Unable to load: %s.%s.%s.%s, exec. time: %d [seconds]\n',...
                knetwks_,kstnms_,kholes_,kcmpnms_,elapsed_time);
            continue;
        end

        %%
        S1orig = S1;
        S1orig = resampleWaveforms(S1orig,finalFs);
        [S1,modelTypes] = transferWaveforms(S1orig,lfc,hfc,npoles,finalFs,units,1,waFlag);

        %%
        if ~isfinite(S1.stla) || ~isfinite(S1.stlo)
            n = n-1;
            elapsed_time = toc;
            fprintf(logFileID,'Unable to load, bad response info: %s.%s.%s.%s, exec. time: %d [seconds]\n',...
                knetwks_,kstnms_,kholes_,kcmpnms_,elapsed_time);
            continue;
        end

        %
        tt = taupTime('iasp91',origdepth,'p,P','km',d_(thisi));
        if isempty(tt)
            n = n-1;
            elapsed_time = toc;
            fprintf(logFileID,'Unable to load: %s.%s.%s.%s, exec. time: %d [seconds]\n',...
                knetwks_,kstnms_,kholes_,kcmpnms_,elapsed_time);
            continue;
        end

        %
        wave_traced_travel_time = min(pull(tt,'time'));
        theoreticalArrival = origt + seconds(wave_traced_travel_time);
        cutStart = theoreticalArrival - minutes(noiseMinutes);
        S1.user0 = wave_traced_travel_time;          % travel time using tables
        S1.eqmag = origmag;
        elapsed_time = toc;
        toc;

        %
        S1 = cutWaveforms(S1,cutStart,0,minutes(durMinutes));
        if isnat(S1.ref)
            fprintf('%s.%s.%s.%s: file is empty, continuing\n',...
                knetwks_,kstnms_,kholes_,kcmpnms_);
            continue;
        end
        if ~strcmpi(units,"WA")
            S1 = scaleWaveforms(detrendWaveforms(S1),100);
        end

        dataVec = S1.d;
        Ntot = length(dataVec);
        duration = seconds(S1.e);

        % write file
        fName_ = sprintf("%s.%s.%s.%s_%s.txt",knetwks_,kstnms_,kholes_,kcmpnms_,datestr(cutStart,30));
        fNameSac_ = sprintf("%s.%s.%s.%s_%s.SAC",knetwks_,kstnms_,kholes_,kcmpnms_,datestr(cutStart,30));
        fNameSAC = fullfile(waveformsDir,fNameSac_);

        %writableChans = ["HN";"BL";"EN";"HH";"SH";"BH"]; %donoterase
        % if strcmp(kcmpnmBase(i),"HN") || strcmp(kcmpnmBase(i),"BL") || ... %donoterase
        %         strcmp(kcmpnmBase(i),"EN") %|| strcmp(kcmpnmstmp(i),"HH") || strcmp(kcmpnmstmp(i),"BH") %donoterase
        try
            disp(fNameSac_);
            sacwrite(fNameSAC,S1);
        catch
            fprintf(2,'couldnt write sac: %s\n',fNameSAC);
        end

        fName = fullfile(waveformsDir,fName_);
        fNameID = fopen(fName,'w');

        %
        if headerVersion == 1
            fuente = "IGEPN";
            modelType = modelTypes(1);
            tipodesuelo = "No Determinado";
            fprintf(fNameID,'#%s\r',"Red Nacional de Acelerógrafos (RENAC) Ecuador - IG-EPN");
            fprintf(fNameID,'#%s\r',"Ladrón de Guevara E11-253, Aptdo. 2759 Quito - Ecuador");
            fprintf(fNameID,'#%s\r',"Teléfonos: (593-2)2225655 ; (593-2)2225627 Fax: (593-2)2567847");
            fprintf(fNameID,'#\r');
            fprintf(fNameID,'#%s\r',"Datos de la estación");
            fprintf(fNameID,'#Nombre de la estación: %s\r',kstnms_);
            fprintf(fNameID,'#Coordenadas de la estación: %3.2f Lat, %3.2f Lon\r',S1.stla,S1.stlo);
            fprintf(fNameID,'#Tipo de suelo: %s\r',tipodesuelo);
            fprintf(fNameID,'#\r');
            fprintf(fNameID,'#%s\r',"Datos del acelerógrafo");
            fprintf(fNameID,'#Modelo del acelerógrafo: %s\r',modelType);
            fprintf(fNameID,'#Frecuencia de muestreo: %d\r',finalFs);
            fprintf(fNameID,'#\r');
            fprintf(fNameID,'#%s\r',"Datos del Sismo");
            fprintf(fNameID,'#Código Evento: %s\r',eventID);
            fprintf(fNameID,'#Fecha/Hora del Sismo (yyyy-mm-dd HH:MM:SS): %s\r',datestr(origt,'yyyy-mm-dd HH:MM:SS.FFF'));
            fprintf(fNameID,'#Magnitud del Sismo: %2.1f\r',origmag);
            fprintf(fNameID,'#Coordenadas del Sismo: %4.3f Lat, %4.3f Lon\r',origlat,origlon);
            fprintf(fNameID,'#Profundidad del Sismo: %4.1f km\r',origdepth);
            fprintf(fNameID,'#Fuente de los datos epicentrales: %s\r',fuente);
            fprintf(fNameID,'#\r');
            fprintf(fNameID,'#%s\r',"Datos de este Registro");
            fprintf(fNameID,'#Componente: %s\r',kcmpnms_);
            fprintf(fNameID,'#Hora de la primera muestra: %s (UTC)\r',datestr(cutStart,'yyyy-mm-dd HH:MM:SS.FFF'));
            fprintf(fNameID,'#Exactitud del tiempo: %4.3f [seg.]\r',1/finalFs);
            fprintf(fNameID,'#Duración del registro: %g [seg.]\r',duration);
            fprintf(fNameID,'#Número total de muestras: %d\r',Ntot);
            fprintf(fNameID,'#Unidades: cm/s^2\r');
            fprintf(fNameID,'#Tipo de filtro: bandpass butterworth, %d poles, lfc: %g [Hz.], hfc: %g [Hz.]\r',npoles,lfc,hfc);
            fprintf(fNameID,'%s\r',"#");
            fprintf(fNameID,'%s\r',"#####################################################################################");
        elseif headerVersion == 2
            % Evento: 201604162358
            %
            % Fecha del evento UTM (aammdd):                  2016    4   16
            %
            % Hora del registro UTM (hhmmss): 23   59 16.00
            %
            % EstaciÛn:                                      AAM2]
            % Componente:                                        E
            %
            % Frecuencia de muestreo (Hz):            1.000000e+02
            %
            % Unidades: cm/s^2
            %
            % _____________________________________________
            fprintf(fNameID,'#%s\n',header1);
            fprintf(fNameID,'#Codigo Evento: %s\n',eventID);
            fprintf(fNameID,'#Tiempo de referencia del archivo: %s\n',datestr(cutStart,'yyyy-mm-dd HH:MM:SS.FFF'));
            fprintf(fNameID,'#Estacion: %s\n',kstnms_);
            fprintf(fNameID,'#Componente: %s\n',kcmpnms_);
            fprintf(fNameID,'#Frecuencia de muestreo (Hz): %d\n',finalFs);
            fprintf(fNameID,'#Latitud, longitud de la Estacion: %3.2f, %3.2f\n',S1.stla,S1.stlo);
            fprintf(fNameID,'#Tiempo de Origen del Evento: %s\n',datestr(origt,'yyyy-mm-dd HH:MM:SS.FFF'));
            fprintf(fNameID,'#Latitud, longitud del Evento: %3.2f, %3.2f\n',origlat,origlon);
            fprintf(fNameID,'#Magnitud del Evento: %2.1f\n',origmag);
            fprintf(fNameID,'#Unidades: cm/s^2\n');
            if ~isfinite(hfc)
                fprintf(fNameID,'#Tipo de filtro: highpass butterworth, %d poles, lfc: %g\n',npoles,lfc);
            else
                fprintf(fNameID,'#Tipo de filtro: bandpass butterworth, %d poles, lfc: %g, hfc: %g\n',npoles,lfc,hfc);
            end
            fprintf(fNameID,'#\n');
        end
        %formatSpec = '%15.8f %15.8f %15.8f %15.8f %15.8f';
        formatSpec = '%f %f %f %f %f';
        lcolumns = 5;
        ldataVec = length(dataVec);
        lrows = ceil(ldataVec/lcolumns);
        ldataVec2 = lrows*lcolumns;
        dataVec = [dataVec; NaN(ldataVec2 - ldataVec,1)]; %#ok<AGROW>

        dataVec = reshape(dataVec,lcolumns,lrows)';
        str = string(compose(formatSpec,dataVec));

        fprintf(fNameID,'%s\r',str);
        fclose(fNameID);
        %end

        %
        S(n,1) = S1;
        fprintf(1,'Done writing: %s.%s.%s.%s, trace: %d, exec. time: %d [seconds]\n',...
            knetwks_,kstnms_,kholes_,kcmpnms_,n,elapsed_time);
        S3(j) = S1;
        k = k + 1;
    end

    %
    if ~k
        disp('i am here');
        continue;
    end

    %
    individualSNCLsProcessed = individualSNCLsProcessed + 1;
    maxAmp = 0;
    kcmpnms__ = char(kcmpnms(locb(1)));
    maxCmp = string(kcmpnms__(end));

    refs = pull(S3,'ref');
    badI = isnat(refs);
    S3 = syncWaveforms(S3(~badI));
    medAmps = median(max(abs(pull(S3))),'omitnan');

    uniqueKstnm = unique(pull(S3,'kstnm'));
    kstnm2(individualSNCLsProcessed,1) = uniqueKstnm(1);

    ncomps = k;
    maxNpts = max(pull(S3,'npts'));
    dVecMagnitude = zeros(maxNpts,1);

    for j = 1:ncomps
        if isnat(S3(j).ref)
            continue;
        end
        dVec = S3(j).d;
        nd = length(dVec);
        dVecMagnitude(1:nd) = dVecMagnitude(1:nd) + dVec.^2;
    end

    dVecMagnitude = sqrt(dVecMagnitude);
    [peakVecMag,peakMagI] = max(abs(dVecMagnitude));
    dVecMagnitudeOrig = dVecMagnitude;
    dVecMagnitude = zpkFilter(dVecMagnitude,-inf,1/2,finalFs,1,true);
    S3(ncomps+1) = S3(ncomps);
    refs = unique(pull(S3,'ref'));
    S3(ncomps+1).ref = refs(1);
    S3(ncomps+1).d = dVecMagnitudeOrig;
    S3(ncomps+1).kstnm = "Envelope";
    S3(ncomps+1).knetwk = "";
    S3(ncomps+1).kcmpnm = "";
    S3(ncomps+1).khole = "";

    %
    tt_ = getTimeVec(S3(1));
    close all;
    [fig,ax] = plotWaveforms(S3,[],[],[],[],[],[],[],Visibility);
    fig.Units = 'normalized';
    fig.OuterPosition = [0 0 0.7 1];

    %
    hold(ax(ncomps+1),'on');
    ll = ax(ncomps+1).Children(1);
    ll.Color(4) = 0.5;
    ll = plot(ax(ncomps+1),tt_,dVecMagnitude,'Linewidth',4);
    ll.Color(4) = 0.75;
    ax(ncomps+1).YLabel.String = unitsStr;
    legend(ax(ncomps+1),'Envelope','Smooth Env.');
    plot(ax(ncomps+1),tt_(peakMagI),peakVecMag,'p','MarkerSize',20,'Linewidth',2);
    text(ax(ncomps+1),tt_(peakMagI)+seconds(2),peakVecMag,sprintf("Max Amplitude: %4.2f, time: %s",...
        peakVecMag,datestr(tt_(peakMagI),'HH:MM:SS.FFF')),'FontSize',fontSize);

    %
    maxAmps = NaN(ncomps,1);
    maxTimes = NaT(ncomps,1);
    for j = 1:ncomps
        if isnat(S3(j).ref)
            continue;
        end

        %
        dVec = S3(j).d;
        signD = sign(dVec);
        dVec = abs(dVec);

        %
        [maxAmp_,maxI] = max(dVec);
        maxAmps(j) = maxAmp_;
        maxTimes(j) = tt_(maxI);
        if maxAmp_ > maxAmp
            maxAmp = maxAmp_;
            kcmpnms__ = char(kcmpnms(locb(j)));
            maxCmp = string(kcmpnms__(end));
        end
        ax(j).YLabel.String = unitsStr;
        hold(ax(j),'on');
        plot(ax(j),tt_(maxI),signD(maxI)*maxAmp_,'p','MarkerSize',20,'Linewidth',2);
        text(ax(j),tt_(maxI)+seconds(2),signD(maxI)*maxAmp_,sprintf("Max Amplitude: %4.2f, time: %s",...
            maxAmp_,datestr(tt_(maxI),'HH:MM:SS.FFF')),'FontSize',fontSize);
    end

    %
    ttI1 = find(tt_ >= theoreticalArrival,1);
    [MAXVECAMP,ttI2] = max(dVecMagnitude);

    triggerThresh = MAXVECAMP/2;
    triggerOnI = find(dVecMagnitude(ttI1:end) >= triggerThresh,1) + ttI1;
    tOn = tt_(triggerOnI);
    triggerOffI = find(flipud(dVecMagnitude(ttI2:end)) >= triggerThresh,1);
    triggerOffI = (maxNpts - triggerOffI + 1); % + ttI2;
    triggerOffI = min([maxNpts triggerOffI]);
    tOff = tt_(triggerOffI);

    fprintf(1,'%s %s %s %f\n',datestr(tOn,'HH:MM:SS.FFF'),datestr(tOff,'HH:MM:SS.FFF'),tOff - tOn,triggerThresh);
    fprintf('\n');
    estimatedDuration = max([minDur seconds(tOff - tOn)]);
    ax(1).Title.String = ...
        {sprintf("Event ID: %s, Origin Lat.: %3.2f, Origin Lon: %3.2f, Magnitude: %2.1f, Depth: %3.2f, Origin Time: %s",...
        eventID,origlat,origlon,origmag,origdepth,datestr(origt,'yyyy-mm-dd HH:MM:SS.FFF')); ...
        sprintf(titleFmt2,...
        sncldist,snclaz,snclbaz,maxAmp,maxCmp,peakVecMag,estimatedDuration)};
    ax(1).Title.FontSize = fontSize+1;

    %
    plotName_ = sprintf("%03d_%s.%s.%s.%s_%s.png",i-1,...
        knetwks_,kstnms_,kholes_,kcmpnms_,datestr(cutStart,30));
    plotName = fullfile(imageDir,plotName_);
    print('-dpng',plotName);

    % print summary information
    for j = 1:ncomps
        if isnat(S3(j).ref)
            continue;
        end
        fprintf(summaryFileID,'%s.%s.%s.%s %g %g %s %g %g %g %g %g\n',...
            knetwks_,kstnms_,kholes_,kcmpnms(locb(j)),...
            maxAmps(j),medAmps,datestr(maxTimes(j),'HH:MM:SS.FFF'),...
            sncldist,snclaz,snclbaz,peakVecMag,estimatedDuration);
    end
    peakAmps(i) = peakVecMag;
    estDurs(i) = estimatedDuration;
end
fclose(logFileID);
fclose(summaryFileID);

%%
S = S(1:n);

close all;
fig = figure('units','normalized','outerposition',[0 0 1/2 1]);
ax = gca;
fig.Visible = Visibility;
hold(ax,'on');

[~,locb] = ismember(uniqSNCLs,allSNCLs);
hypocentralDistances = d_(locb);

badI = ~isfinite(estDurs) | ~isfinite(peakAmps) | ~isfinite(hypocentralDistances);
%if sum(badI) >= individualSNCLsProcessed %<-- i dont know why i wrote this option
if sum(~badI) < 2
    %fprintf(2,"something went wrong, too many bad events.\n");
    fprintf(2,"something went wrong, not enough good events.\n");
    return;
end

kstnm2 = kstnm2(1:individualSNCLsProcessed);
estDurs(badI) = [];
peakAmps(badI) = [];
hypocentralDistances(badI) = [];

%d_ = sqrt(d_.^2 + origdepth^2);
%disp(d_);
minDist = 10.^floor(log10(min(hypocentralDistances)));
maxDist = 10.^ceil(log10(max(hypocentralDistances)));

%
plot(ax,hypocentralDistances,peakAmps,'.','MarkerSize',25);
xlim([minDist maxDist]);

grid on;
xlabel('Hypocentral Distance [km.]');
ylabel(strcat('Peak Vector Amplitude, [',unitsStr,']'));
if individualSNCLsProcessed > 2
    b_ = flipud(robustfit(log10(hypocentralDistances),log10(peakAmps)));
    yq = polyval(b_,log10(hypocentralDistances));
    ll = plot(ax,hypocentralDistances,10.^yq,'linewidth',3);
    ll.Color(4) = 0.5;
end

%
title(ax,['slope: ',num2str(b_(1))]);
text(ax,hypocentralDistances,peakAmps,kstnm2,'FontSize',12);
plotName_ = sprintf("%s_%s.png","AmplitudeVsDistance",upper(units));
plotName = fullfile(writeDir,plotName_);
ax.XScale = 'log';
ax.YScale = 'log';
print('-dpng',plotName);

close all;
fig = figure('units','normalized','outerposition',[0 0 1/2 1]);
ax = gca;
fig.Visible = Visibility;
hold(ax,'on');

plot(ax,hypocentralDistances,estDurs,'.','MarkerSize',25);
xlim([minDist maxDist]);

grid on;
xlabel('Hypocentral Distance [km.]');
ylabel('Estimated Duration [sec.]');
if individualSNCLsProcessed > 2
    b_ = flipud(robustfit(log10(hypocentralDistances),log10(estDurs)));
    yq = polyval(b_,log10(hypocentralDistances));
    ll = plot(ax,hypocentralDistances,10.^yq,'linewidth',3);
    ll.Color(4) = 0.5;
end

%
title(ax,['slope: ',num2str(b_(1))]);
text(ax,hypocentralDistances,estDurs,kstnm2,'FontSize',12);
plotName_ = sprintf("%s_%s.png","DurationVsDistance",upper(units));
plotName = fullfile(writeDir,plotName_);
ax.XScale = 'log';
ax.YScale = 'log';

print('-dpng',plotName);
pause(2);
close all;