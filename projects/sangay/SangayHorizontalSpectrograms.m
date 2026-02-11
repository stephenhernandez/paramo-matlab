function S = SangayHorizontalSpectrograms()
cd ~/research/now/sangay/horizontalSpectrograms/

npoles = 2;
newFs = 40;

nDays = 2;
kstnm = "SAGA";
kcmpnms = ["HHZ";"HHN";"HHE";"BDF"];

%%
lfc = 1/1000;
hfc = -inf;
rotAng = 35;
units = 'acc';
tw = 0.001;
value = 0;

%%
dayVec = (datetime(2022,10,20):-1:datetime(2022,09,24))';
lDays = length(dayVec);

%%
tic;
for i = 1:lDays
    tic;
    dayStart = dayVec(i);
    fprintf('attempting to process: %s\n',datestr(dayStart));
    S = loadWaveforms(dayStart,nDays,kstnm,kcmpnms,"EC",["";"01"]);

    %
    refs = pull(S,'ref');
    badI = isnat(refs);
    S(badI) = [];
    lS = length(S);

    if lS < 4
        fprintf('skipping: %s\n',datestr(dayStart));
        continue;
    end

    if isnat(S(1).ref)
        fprintf('skipping: %s\n',datestr(dayStart));
        continue;
    end
    
    %%
    fprintf('processing: %s\n',datestr(dayStart));

    %%
    S = intWaveforms(detrendWaveforms(differentiateWaveforms(S)));
    S = taperWaveforms(S,tw);
    S = transferWaveforms(S,lfc,hfc,npoles,100,units,true);
    S(1:3) = scaleWaveforms(S(1:3),1e9);
    S(4) = scaleWaveforms(S(4),1e6);
    S = rotateWaveforms(S,rotAng,3,2,true);
    S = resampleWaveforms(S,newFs);
    S = syncWaveforms(S);
    S = interpolateWaveforms(S,value);

    %%
    for j = 1:lS
        Sj = S(j);
        if isnat(Sj.ref)
            continue;
        end
        
        charkcmpnm = char(Sj.kcmpnm);
        if ~exist(charkcmpnm,'dir')
            mkdir(charkcmpnm)
        end
        cd(charkcmpnm);
        
        if ~exist('VLF-HF','dir')
            mkdir('VLF-HF')
        end
        cd('VLF-HF');
        
        %%
        logSpectrogram(Sj,512,16,3/4,1/500,20);
        clim([20 180]);
        if j == 4
            clim([40 180]);
        end
        
        colormap parula;
        fName = ['EC.SAGA..',charkcmpnm,'_',datestr(dayStart,'yyyy.mm.dd'),'_VLF-HF.jpg'];
        print('-djpeg',fName);
        cd ..;
        close all;
        cd ..
        toc;
    end
end
