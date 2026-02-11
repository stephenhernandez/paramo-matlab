function [G,d] = tensorInversion6(O,t0,testLat,testLon,testDepth,th,lfc,hfc,npoles,secDur,units)
lS = length(O);
O = cutSacData(O,t0,0,secDur);
refEllipse = referenceEllipsoid('wgs84');
stla = pull(O,'stla');
stlo = pull(O,'stlo');

%%
G = [];
d = G;
[dist,azs] = distance(testLat,testLon,stla,stlo,refEllipse);
[~,bazs] = distance(stla,stlo,testLat,testLon,refEllipse);
dist = dist*1e-3;

%% get synthetics
for i = 1:lS
    d = [d; pull(O(i))];     %#ok<AGROW>
    
    %%
    Fs = round(1./O(i).delta);
    cmp = char(O(i).kcmpnm);
    cmp = cmp(end);
    disp(dist(i))
    [T,R,Z] = getGreensFunctionsAll(round(dist(i)),round(testDepth),azs(i)); %,1e0,'~/gf/nicoyaGF/');
    stf = getSTF(th,Fs);
    
    if strcmp(cmp,'Z')
        Z = filterWaveforms(Z,lfc,hfc,npoles);
        
        %%
        if strcmp(units,'vel')
            Z = differentiateWaveforms(Z);
        end
        
        if strcmp(units,'acc')
            Z = differentiateWaveforms(Z,2);
        end
        
        Z = cutWaveforms(Z,datetime(1970,01,01),0,secDur);
        Z = resampleWaveforms(Z,Fs);
        Z = convWaveforms(Z,stf);
        Z = pull(Z,'d');
        
        G = [G; Z];                         %#ok<AGROW>
    else
        %% filter
        T = filterWaveforms(T,lfc,hfc,npoles);
        R = filterWaveforms(R,lfc,hfc,npoles);
        
        %%
        if strcmp(units,'vel')
            T = differentiateWaveforms(T);
            R = differentiateWaveforms(R);
        end
        
        if strcmp(units,'acc')
            T = differentiateWaveforms(T,2);
            R = differentiateWaveforms(R,2);
        end
        
        %% cut
        T = cutWaveforms(T,datetime(1970,01,01),0,secDur);
        R = cutWaveforms(R,datetime(1970,01,01),0,secDur);
        
        %% resample
        T = resampleWaveforms(T,Fs);
        R = resampleWaveforms(R,Fs);
        
        %% convolve
        T = convWaveforms(T,stf);
        R = convWaveforms(R,stf);
        
        %% extract synthetics
        T = pull(T,'d');
        R = pull(R,'d');
        
        %%
        baz = bazs(i);
        if baz >= 360
            baz = baz - 360;
        end
        
        if baz <= 180
            rotang = baz + 180;
        else
            rotang = baz - 180;
        end
        
        %%
        for ii = 1:5
            ttemp = T(:,ii);
            rtemp = R(:,ii);
            [ee,nn] = rotate2d(ttemp,rtemp,rotang);
            R(:,ii) = nn;
            T(:,ii) = ee;
        end
        
        %%
        if strcmp(cmp,'E')
            G = [G; T];                     %#ok<AGROW>
        else
            G = [G; R];                     %#ok<AGROW>
        end
    end
end

%% perform inversion
%m = ((G'*G)^-1)*G'*d;
%synth = G*m;
%m = [m(1) m(2) -(m(1)+m(2)) m(3) m(4) m(5)];
%mw = mt2mw(m);
