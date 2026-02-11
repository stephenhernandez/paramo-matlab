function [MLv,nMLv,mB,nmB,ML,nML,Mjma,nMjma,MsBB,nMsBB,Mwp,nMwp,mb,nmb,magerr2] = readMagnitudeInformation(event_file)
fid = fopen(event_file);
eofstat = feof(fid);

%preallocate
% MLv = populateMagnitudeStruct(10);
% ML = MLv;
% Mjma = MLv;
% MsBB = MLv;
% Mwp = MLv;
% mB = MLv;
% mb = MLv;

%%
sCount = 0;
mFlag = false;
nMLv = 0;
nML = 0;
nMjma = 0;
nMsBB = 0;
nMwp = 0;
nmB = 0;
nmb = 0;
magerr2 = NaN;

%get the job done
while ~eofstat
    buff = fgetl(fid);
    if buff ~= -1
        tmp = sscanf(buff,'%d %s');
        charstr = char(tmp(2:end))';
        if strcmp(charstr,'Phase')
            mFlag = false;
        elseif strcmp(charstr,'Station')
            nmag = tmp(1);
            MLv = populateMagnitudeStruct(nmag);
            ML = MLv;
            Mjma = MLv;
            MsBB = MLv;
            Mwp = MLv;
            mB = MLv;
            mb = MLv;
            mFlag = true;
        elseif mFlag
            sCount = sCount + 1;
            if sCount > 0 && sCount < nmag + 1
                C = textscan(buff,'%s %s %s %f %f %s %f %f %f');
                type = C{6};
                if strcmp(type,'MLv')
                    nMLv = nMLv+1;
                    stnm = C{1}; MLv(nMLv).stnm = string(stnm{1});
                    ntwk = C{2}; MLv(nMLv).ntwk = string(ntwk{1});
                    chan = C{3}; chan = chan{1}; MLv(nMLv).chan = string(chan(1:end));
                    MLv(nMLv).dist = C{4};
                    MLv(nMLv).azimuth = C{5};
                    MLv(nMLv).type = string(type);
                    MLv(nMLv).value = C{7};
                    MLv(nMLv).res = C{8};
                    amp_ = C{9};
                    if isempty(amp_)
                        MLv(nMLv).amp = NaN;
                    else
                        MLv(nMLv).amp = amp_;
                    end
                    
                elseif strcmp(type,'mB')
                    nmB = nmB+1;
                    stnm = C{1}; mB(nmB).stnm = string(stnm{1});
                    ntwk = C{2}; mB(nmB).ntwk = string(ntwk{1});
                    chan = C{3}; chan = chan{1}; mB(nmB).chan = string(chan(1:end));
                    mB(nmB).dist = C{4};
                    mB(nmB).azimuth = C{5};
                    mB(nmB).type = string(type);
                    mB(nmB).value = C{7};
                    mB(nmB).res = C{8};
                    mB(nmB).amp = C{9};
                    
                elseif strcmp(type,'mb')
                    nmb = nmb+1;
                    stnm = C{1}; mb(nmb).stnm = string(stnm{1});
                    ntwk = C{2}; mb(nmb).ntwk = string(ntwk{1});
                    chan = C{3}; chan = chan{1}; mb(nmb).chan = string(chan(1:end));
                    mb(nmb).dist = C{4};
                    mb(nmb).azimuth = C{5};
                    mb(nmb).type = string(type);
                    mb(nmb).value = C{7};
                    mb(nmb).res = C{8};
                    mb(nmb).amp = C{9};
                    
                elseif strcmp(type,'Mjma')
                    nMjma = nMjma+1;
                    stnm = C{1}; Mjma(nMjma).stnm = string(stnm{1});
                    ntwk = C{2}; Mjma(nMjma).ntwk = string(ntwk{1});
                    chan = C{3}; chan = chan{1}; Mjma(nMjma).chan = string(chan(1:end));
                    Mjma(nMjma).dist = C{4};
                    Mjma(nMjma).azimuth = C{5};
                    Mjma(nMjma).type = string(type);
                    Mjma(nMjma).value = C{7};
                    Mjma(nMjma).res = C{8};
                    Mjma(nMjma).amp = C{9};
                    Mjma = Mjma';
                    
                elseif strcmp(type,'ML')
                    nML = nML+1;
                    stnm = C{1}; ML(nML).stnm = string(stnm{1});
                    ntwk = C{2}; ML(nML).ntwk = string(ntwk{1});
                    chan = C{3}; chan = chan{1}; ML(nML).chan = string(chan(1:end));
                    ML(nML).dist = C{4};
                    ML(nML).azimuth = C{5};
                    ML(nML).type = string(type);
                    ML(nML).value = C{7};
                    ML(nML).res = C{8};
                    ML(nML).amp = C{9};
                    
                elseif strcmp(type,'Ms(BB)')
                    nMsBB = nMsBB+1;
                    stnm = C{1}; MsBB(nMsBB).stnm = string(stnm{1});
                    ntwk = C{2}; MsBB(nMsBB).ntwk = string(ntwk{1});
                    chan = C{3}; chan = chan{1}; MsBB(nMsBB).chan = string(chan(1:end));
                    MsBB(nMsBB).dist = C{4};
                    MsBB(nMsBB).azimuth = C{5};
                    MsBB(nMsBB).type = string(type);
                    MsBB(nMsBB).value = C{7};
                    MsBB(nMsBB).res = C{8};
                    MsBB(nMsBB).amp = C{9};
                    
                elseif strcmp(type,'Mwp')
                    nMwp = nMwp+1;
                    stnm = C{1}; Mwp(nMwp).stnm = string(stnm{1});
                    ntwk = C{2}; Mwp(nMwp).ntwk = string(ntwk{1});
                    chan = C{3}; chan = chan{1}; Mwp(nMwp).chan = string(chan(1:end));
                    Mwp(nMwp).dist = C{4};
                    Mwp(nMwp).azimuth = C{5};
                    Mwp(nMwp).type = string(type);
                    Mwp(nMwp).value = C{7};
                    Mwp(nMwp).res = C{8};
                    Mwp(nMwp).amp = C{9};
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
mB = mB(1:nmB);
mb = mb(1:nmb);
if nMLv > 2
    MLv_values = pull(MLv,'value');
    magerr2 = mad(MLv_values,1);
end
