function [S,f,fI] = fdWhiten(S,lfc,hfc,dW,Fs,timeDomain,verboseFlag,plotFlag)
if nargin < 2; lfc = -inf; end
if nargin < 3; hfc = -inf; end
if nargin < 4; dW = 0; end
if nargin < 6; timeDomain = true; end
if nargin < 7; verboseFlag = true; end
if nargin < 8; plotFlag = false; end

%% frequency domain normalization
isstructS = isstruct(S);
cornersfin = isfinite([lfc hfc]);
npoles = 1;
zeroPhaseFlag = true;

%%
if isstructS
    lS = length(S);
    for i = 1:lS
        n = S(i).npts;
        Fs = 1./S(i).delta;
        d = S(i).d;
        D = fft(d);
        df = Fs/n;

        %%
        fI = true(n,1);
        if any(cornersfin)
            f = df*(0:n-1)';
            if cornersfin(1)    % high-pass filter
                fI(f < lfc | f > Fs-lfc) = false;
            end

            if cornersfin(2)    % low-pass filter
                fI(f > hfc & f < Fs-hfc) = false;
            end
        end

        %%
        normWeight = abs(D);
        if dW > 0
            N = ceil(dW/df);
            if N > 3
                if verboseFlag
                    fprintf("smoothing the spectrum with a %d-point window\n",N);
                end
                normWeight = medfiltSH(normWeight,N,zeroPhaseFlag);
                %normWeight = zpkFilter(normWeight,-inf,1/N,1,npoles,zeroPhaseFlag);
                %figure(); semilogy([normWeight normWeight2],'.'); zoom on; grid on;
            end
        end
        fI = fI & normWeight > 0;
        D = D./normWeight;
        D(~fI) = 0;


        %%
        d = ifft(D,'symmetric');
        drms = rms(d);
        S(i).d = d/drms;

        %%
        [minVals,maxVals,meanVals] = minmaxmean(d);
        S(i).depmin = minVals;
        S(i).depmax = maxVals;
        S(i).depmen = meanVals;
    end
else
    if nargin < 5; Fs = 100; end
    
    %%
    % input data are not waveform structs. treat as matrix where columns
    % are individual traces.
    n = size(S,1);
    if timeDomain
        % convert to frequency domain
        S = fft(S,n);
    end
    df = Fs/n;
    
    %%
    f = df*(0:n-1)';
    fI = true(n,1);
    if any(cornersfin)
        if cornersfin(1)
            fI(f < lfc | f > Fs-lfc) = false;
            %fI(f > hfc) = false;
        end

        if cornersfin(2)
            fI(f > hfc & f < Fs-hfc) = false;
            %fI(f < lfc) = false;
        end
    end
    
    %%
    normWeight = abs(S);
    if dW > 0
        N = ceil(dW/df);
        if verboseFlag
            fprintf('smoothing the spectrum with a %d-point window: %d\n',N,n);
        end

        convSmoothingFlag = false;
        medfiltSmoothingFlag = true;
        if N > 1
            if convSmoothingFlag
                normWeight = fftfilt(ones(N,1)/N,normWeight);
                normWeight = flipud(normWeight);
                normWeight = fftfilt(ones(N,1)/N,normWeight);
                normWeight = flipud(normWeight);
            elseif medfiltSmoothingFlag
                normWeight = medfiltSH(normWeight,N,zeroPhaseFlag);
            else
                fprintf('not applying smoothing, something went wrong\n');
            end
        end
    end
    fI = fI & normWeight > 0;

    if plotFlag
        Sorig = S;
    end
    S = S./normWeight;
    S(~fI) = 0;

    if plotFlag
        figure();
        ax_(1) = subplot(311);
        semilogx(f,abs(Sorig(:,1))); zoom on; grid on; 
        hold on;
        semilogx(f,normWeight(:,1));
        ax_(2) = subplot(312);
        semilogx(f,abs(Sorig(:,1)./normWeight(:,1)));
        zoom on; grid on;
        ax_(3) = subplot(313);
        semilogx(f,abs(S));
        zoom on; grid on;
        linkaxes(ax_,'x');
    end

    if timeDomain
        % convert back to time domain
        S = ifft(S,[],1,'symmetric');
    end
end
