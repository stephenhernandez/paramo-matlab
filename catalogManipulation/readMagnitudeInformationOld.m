function [origmag,nMLv,MLv,nML,ML,...
    nMjma,Mjma,nMsBB,MsBB,nMwp,Mwp] = readMagnitudeInformationOld(event_file)
fid = fopen(event_file);
eofstat = feof(fid);

%preallocate
origmag = NaN;
magCount = 0;
sCount = 0;
nFlag = false;
mFlag = false;
nMLv = 0;
nML = 0;
nMjma = 0;
nMsBB = 0;
nMwp = 0;

%get the job done
while ~eofstat
    buff = fgetl(fid);
    if buff ~= -1
        tmp = sscanf(buff,'%d %s');
        charstr = char(tmp(2:end))';
        if strcmp(charstr,'Network')
            nmags = tmp(1);
            nFlag = true;
        elseif strcmp(charstr,'Phase')
            nFlag = false;
        elseif strcmp(charstr,'Station')
            nmag = tmp(1);
            MLv = populateMagnitudeStruct(nmag);
            ML = MLv;
            Mjma = MLv;
            MsBB = MLv;
            Mwp = MLv;
            mFlag = true;
        elseif nFlag
            magCount = magCount + 1;
            if magCount <= nmags
                C = textscan(buff,'%s %f %s %s %s %s');
                if strcmp(C{4},'preferred')
                    %disp('4th column preferred')
                    origmag = C{2};
                elseif strcmp(C{5},'preferred')
                    %disp('5th column preferred')
                    origmag = C{2};
                elseif strcmp(C{6},'preferred')
                    %disp('6th column preferred')
                    origmag = C{2};
                end
            end
        elseif mFlag
            sCount = sCount + 1;
            if sCount > 1 && sCount < nmag + 2
                C = textscan(buff,'%s %s %f %d %s %f %f %f');
                type = C{5};
                if strcmp(type,'MLv')
                    nMLv = nMLv+1;
                    %disp(['MLv: ',num2str(nMLv)]);
                    stnm = C{1}; MLv(nMLv).stnm = string(stnm{1});
                    ntwk = C{2}; MLv(nMLv).ntwk = string(ntwk{1});
                    MLv(nMLv).dist = C{3};
                    MLv(nMLv).azimuth = C{4};
                    MLv(nMLv).type = string(type);
                    MLv(nMLv).value = C{6};
                    MLv(nMLv).res = C{7};
                    MLv(nMLv).amp = C{8};
                elseif strcmp(type,'Mjma')
                    nMjma = nMjma+1;
                    %disp(['Mjma: ',num2str(nMjma)]);
                    stnm = C{1}; Mjma(nMjma).stnm = string(stnm{1});
                    ntwk = C{2}; Mjma(nMjma).ntwk = string(ntwk{1});
                    Mjma(nMjma).dist = C{3};
                    Mjma(nMjma).azimuth = C{4};
                    Mjma(nMjma).type = string(type);
                    Mjma(nMjma).value = C{6};
                    Mjma(nMjma).res = C{7};
                    Mjma(nMjma).amp = C{8};
                elseif strcmp(type,'ML')
                    nML = nML+1;
                    %disp(['ML: ',num2str(nML)]);
                    stnm = C{1}; ML(nML).stnm = string(stnm{1});
                    ntwk = C{2}; ML(nML).ntwk = string(ntwk{1});
                    ML(nML).dist = C{3};
                    ML(nML).azimuth = C{4};
                    ML(nML).type = string(type);
                    ML(nML).value = C{6};
                    ML(nML).res = C{7};
                    ML(nML).amp = C{8};
                elseif strcmp(type,'Ms(BB)')
                    nMsBB = nMsBB+1;
                    %disp(['MsBB: ',num2str(nMsBB)]);
                    stnm = C{1}; MsBB(nMsBB).stnm = string(stnm{1});
                    ntwk = C{2}; MsBB(nMsBB).ntwk = string(ntwk{1});
                    MsBB(nMsBB).dist = C{3};
                    MsBB(nMsBB).azimuth = C{4};
                    MsBB(nMsBB).type = string(type);
                    MsBB(nMsBB).value = C{6};
                    MsBB(nMsBB).res = C{7};
                    MsBB(nMsBB).amp = C{8};
                elseif strcmp(type,'Mwp')
                    nMwp = nMwp+1;
                    %disp(['Mwp: ',num2str(nMwp)])
                    stnm = C{1}; Mwp(nMwp).stnm = string(stnm{1});
                    ntwk = C{2}; Mwp(nMwp).ntwk = string(ntwk{1});
                    Mwp(nMwp).dist = C{3};
                    Mwp(nMwp).azimuth = C{4};
                    Mwp(nMwp).type = string(type);
                    Mwp(nMwp).value = C{6};
                    Mwp(nMwp).res = C{7};
                    Mwp(nMwp).amp = C{8};
                end
            end
        end
    end
    eofstat = feof(fid);
end
fclose(fid);

%%
MLv = MLv(1:nMLv);
ML = ML(1:nML);
Mjma = Mjma(1:nMjma);
MsBB = MsBB(1:nMsBB);
Mwp = Mwp(1:nMwp);