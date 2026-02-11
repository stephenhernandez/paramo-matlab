%function regionalReventadorSensorsForPablo()
kstnm = ["REVS";"CAYR";"CASC";"ANTI";"ANTS"];
%kcmpnm = ["HHZ";"SHZ";"HHZ";"HHZ";"HHZ"];

%%
%cd ~/observatorios/07-REVENTADOR/datos/datos-regionales
tStart = datetime(2020,09,01);
tEnd = datetime(2020,09,09);
days = tStart:tEnd;
lDays = length(days)-1;

%%
lfc = 0.5;
hfc = 2;

%%
for i = 1:lDays
    S = loadWaveforms(days(i),1,kstnm,["HHZ";"SHZ"],"EC");
    %cd ~/observatorios/07-REVENTADOR/datos/datos-regionales;
    cd ~/
    tStart_ = dateshift(days(i),'start','day') + hours(0) + minutes(0);
    littleWindows = (tStart_:minutes(20):tStart_+1)';
    ll = length(littleWindows)-1;
    if ~isempty(S)
        lS = length(S);
        [~,mm,~,~,~,~] = datevec(tStart_);
        if mm < 10
            cd(['2020',num2str(mm)]);
        else
            disp('script is broken');
            break;
        end
        
        %%
        for j = 1:lS
            S_ = filterWaveforms(detrendWaveforms(S(j)),lfc,hfc);
            for k = 1:ll
                littleWindows_ = littleWindows(k);
                Scut = cutWaveforms(S_,littleWindows_,0,minutes(20),true);
                [~,~,dd,HH,MM,SS] = datevec(littleWindows_);
                cdDir = getdayhourminsec(dd,HH,MM,SS);
                if ~isnat(Scut.ref)
                    %%
                    disp(cdDir);
                    if exist(cdDir,'dir')
                        cd(cdDir);
                    else
                        mkdir(cdDir)
                        ls
                        cd(cdDir);
                    end
                    
                    %%
                    kstnm_ = Scut.kstnm;
                    chan_ = Scut.kcmpnm;
                    kntwk = "EC";
                    fileName = strcat(char(kstnm_),'.',char(chan_),'.',char(kntwk),'.sac');
                    disp(fileName)
                    sacwrite(fileName,Scut);
                    cd ..;
                end
            end
        end
        cd ..;
    end
end


%%
function cdDir = getdayhourminsec(dd,HH,MM,SS)
if dd < 10
    cdDir = ['0',num2str(dd)];
    if HH < 10
        cdDir = [cdDir,'0',num2str(HH)];
        if MM < 10
            cdDir = [cdDir,'0',num2str(MM)];
            if SS < 10
                cdDir = [cdDir,'0',num2str(SS)];
            else
                cdDir = [cdDir,num2str(SS)];
            end
        else
            cdDir = [cdDir,num2str(MM)];
            if SS < 10
                cdDir = [cdDir,'0',num2str(SS)];
            else
                cdDir = [cdDir,num2str(SS)];
            end
        end
    else
        cdDir = [cdDir,num2str(HH)];
        if MM < 10
            cdDir = [cdDir,'0',num2str(MM)];
            if SS < 10
                cdDir = [cdDir,'0',num2str(SS)];
            else
                cdDir = [cdDir,num2str(SS)];
            end
        else
            cdDir = [cdDir,num2str(MM)];
            if SS < 10
                cdDir = [cdDir,'0',num2str(SS)];
            else
                cdDir = [cdDir,num2str(SS)];
            end
        end
    end
else
    cdDir = num2str(dd);
    if HH < 10
        cdDir = [cdDir,'0',num2str(HH)];
        if MM < 10
            cdDir = [cdDir,'0',num2str(MM)];
            if SS < 10
                cdDir = [cdDir,'0',num2str(SS)];
            else
                cdDir = [cdDir,num2str(SS)];
            end
        else
            cdDir = [cdDir,num2str(MM)];
            if SS < 10
                cdDir = [cdDir,'0',num2str(SS)];
            else
                cdDir = [cdDir,num2str(SS)];
            end
        end
    else
        cdDir = [cdDir,num2str(HH)];
        if MM < 10
            cdDir = [cdDir,'0',num2str(MM)];
            if SS < 10
                cdDir = [cdDir,'0',num2str(SS)];
            else
                cdDir = [cdDir,num2str(SS)];
            end
        else
            cdDir = [cdDir,num2str(MM)];
            if SS < 10
                cdDir = [cdDir,'0',num2str(SS)];
            else
                cdDir = [cdDir,num2str(SS)];
            end
        end
    end
end
end