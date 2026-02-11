function [uPos,gPos,cc,u,residualRMS,lastmaxcc] = iterative_deconvolution(u,stf,max_nPos,plotFlag)%,thresh) %stfRange)
if nargin < 3
    max_nPos = 300;
end

if nargin < 4
    plotFlag = false;
end

%%
u = detrend(u);
u = u/rssq(u);
winlen = length(stf);
g = zeros(size(u));
uOrig = u;

%% iteration 1
maxIter = 5;

residualRMS = NaN(maxIter,1);
counter = 0;
posCounter = 0;
lastmaxcc = 1;
box = ones(winlen,1);

%%
origVersion = ~true; flag2 = true; %both true == original algorithm but slow
if origVersion
    while counter < maxIter
        %disp([posCounter lastmaxcc])
        counter = counter + 1;
        cc = flipud(fftfilt(stf,flipud(u)));
        normers = sqrt(abs(flipud(fftfilt(box,flipud(u.^2)))));
        normcc = cc./normers;
        %signcc = sign(cc);
        if flag2
            [maxCC,maxI] = max(abs(cc));
            w = normcc(maxI);
            %signcc_ = signcc(maxI);
            g(maxI) = g(maxI) + w.*maxCC;

            lastmaxcc = maxCC;
            u2 = fftfilt(stf,g);

            u = uOrig - u2;
            residualRMS(counter) = rms(cc); %%rms(u);

            if sign(w) < 0
                continue;
            end

            posCounter = posCounter + 1;
            %disp([posCounter counter])
            if posCounter == max_nPos
                counter = maxIter;
            end
        else
            [pks,locb] = findpeaks(abs(cc),'MinPeakDistance',15);
            [~,sortPKS] = sort(pks,'descend');
            locb = locb(sortPKS);
            llocb = length(locb);
            signcc = sign(cc);

            locI = min([llocb max_nPos]);
            maxI = locb(1:locI);
            maxCC = cc(maxI);
            w = abs(normcc(maxI));

            signcc_ = signcc(maxI) > 0;
            g(maxI) = g(maxI) + w.*maxCC;

            u2 = fftfilt(stf,g);
            u = uOrig - u2;
            residualRMS(counter) = rms(cc); %%rms(u);

            if ~any(signcc_)
                continue;
            end

            posCounter = posCounter + 1;
            %disp([posCounter counter])
        end
    end

    %plotFlag = false;
    if plotFlag
        plotter_(cc,g,u,u2,uOrig,stf)
        figure();
        plot(residualRMS,'o'); zoom on; grid on;
    end

    gPos = g;
    gPos(gPos<0) = 0;
    %gPos = sign(gPos); %experimental!
    uPos = fftfilt(stf,gPos);
else
    while posCounter < 1%max_nPos
        counter = counter + 1;
        rawcc = flipud(fftfilt(stf,flipud(u)));
        normers = sqrt(abs(flipud(fftfilt(box,flipud(u.^2)))));
        cc = rawcc./normers;
        %[~,locb] = findpeaks(abs(cc),'MinPeakDistance',15,'MinPeakHeight',0.05);
        %[~,locb] = findpeaks((cc),'MinPeakDistance',15,'MinPeakHeight',0.3);
        [pks,locb] = findpeaks((cc),'MinPeakDistance',32*30,'MinPeakHeight',0.1);
        %rawpks = rawcc(locb);
        [~,sortPKS] = sort(pks,'descend');
        locb = locb(sortPKS);
        llocb = length(locb);
        signcc = sign(cc);

        locI = min([llocb max_nPos]);
        maxI = locb(1:locI);
        maxCC = cc(maxI); %rawcc(maxI);

        signcc_ = signcc(maxI);
        %g(maxI) = g(maxI) + signcc_*rawcc(maxI);
        g(maxI) = g(maxI) + maxCC;

        if any(signcc_) %&& lastmaxcc >= thresh
            posCounter = posCounter + 1;
            disp([posCounter locI])
            lastmaxcc = maxCC;
            u2 = fftfilt(stf,g);

            u = uOrig - u2;
            residualRMS(counter) = rms(u);
        end

        if counter == maxIter
            posCounter = max_nPos;
        end
    end

    %plotFlag = false;
    if plotFlag
        plotter_(cc,g,u,u2,uOrig,stf)
        figure();
        plot(residualRMS,'o'); zoom on; grid on;
    end

    gPos = g;
    gPos(gPos<0) = 0;
    %gPos = sign(gPos); %experimental!
    uPos = fftfilt(stf,gPos);
end
end

function plotter_(cc,g,u,u2,uOrig,stf)
figure('units','normalized','outerposition',[0.1 0 0.8 1]);
ax(1) = subplot(421);
plot(cc); zoom on;
title('cc');

ax(2,1) = subplot(422);
stem(g,'.-'); zoom on;

ax(3,1) = subplot(423);
plot(uOrig,'linewidth',2); hold on;
plot(u2);
title('synthetic');

ax(4,1) = subplot(424);
gPos = g;
gPos(gPos<0) = 0;
synth2 = fftfilt(stf,gPos);
plot(uOrig,'k','linewidth',0.2); hold on;
ax(4,1).ColorOrderIndex = 1;
plot(synth2,'linewidth',2);
title('g-pos')

ax(5,1) = subplot(425);
plot(u);
title('previous residual');

ax(6,1) = subplot(426);
gNeg = g;
gNeg(gNeg>=0) = 0;
synth3 = fftfilt(stf,gNeg);
plot(synth3); %[synth2 synth3]);
title('g-neg');
%legend('g-pos','g-neg')

ax(7,1) = subplot(427);
plot(uOrig-u2);
title('current residual');

ax(8,1) = subplot(428);
%plot(uOrig,'linewidth',2); hold on;
plot(u2,'linewidth',1);
title('synthetic');

linkaxes(ax,'xy');
ylim(0.05*[-1 1])
end
