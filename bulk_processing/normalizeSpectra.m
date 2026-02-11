function [pxx,tcut,h,ax,fxx] = normalizeSpectra(S,secDur,nOverlap,pmtmFlag,plotFlag,minLev,normFlag,maxLev)
if nargin < 2; secDur = 0.5*60; end
if nargin < 3; nOverlap = 0.5; end
if nargin < 4; pmtmFlag = true; end
if nargin < 5; plotFlag = true; end
if nargin < 6; minLev = -100; end
if nargin < 7; normFlag = true; end
if nargin < 8; maxLev = 100; end

lS = size(S,1);
h = gobjects(lS,1);
ax = gobjects(lS,1);
for i = 1:lS
    Scut = S(i);
    Fs = round(1./Scut.delta);
    if nOverlap > 1
        nOverlap = round(nOverlap);
    else
        nOverlap = round(nOverlap*secDur*Fs);
    end
    
    t = getTimeVec(Scut);
    [dcut,~,indices] = cutWindows(Scut.d,secDur*Fs,nOverlap,true);
    tcut = t(indices);
    
    ld = size(dcut,1);
    npow2 = 2+nextpow2(ld);
    nfft = 2^npow2;
    
    if pmtmFlag
        [pxx,fxx] = pmtm(dcut.*parzenwin(ld),1.5,nfft,Fs);
    else
        [pxx,fxx] = pwelch(dcut.*parzenwin(ld),[],[],nfft,Fs);
    end
    
    if normFlag
        pxx = bsxfun(@rdivide,pxx,max(abs(pxx)));
    end
    pxx = 10*log10(pxx.^2);
    
    if plotFlag
        figure('units','normalized','outerposition',[0.2 0.2 1 0.8]);
        ax(i) = gca;
        pxx(pxx <= minLev) = NaN;
        h(i) = imagesc(datenum(tcut),fxx,pxx);
        axis xy;
        colorbar;
        if normFlag
            caxis(ax(i),[minLev 0]);
        else
            caxis(ax(i),[minLev maxLev]);
        end
        datetick('x','keeplimits');
        zoom on;
        set(h(i), 'alphadata', isfinite(pxx));
        title(strcat(Scut.kstnm,", ",Scut.kcmpnm,", ",datestr(Scut.ref,1)));
    end
end