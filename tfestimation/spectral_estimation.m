function [AmpSmooth,CrossSpec,PhaseSmooth] = spectral_estimation(S,Hbu,...
    totN,dpssSeq,nfft,lambda,overlapPrcnt,detrendFlag,...
    GAINCORRECTION,HILBERTFLAG,CANUSEGPU)
%
% [decon,phase,xcohe,mscohe,Sxy,Sxx,Syy] = fd_decon(Sf,Hbu,totN,dpssSeq,nfft,...
% overlapPrcnt,detrendFlag,GAINCORRECTION,HILBERTFLAG,CANUSEGPU)
%

%%
if nargin < 8
    detrendFlag = true;
end

if nargin < 9
    GAINCORRECTION = false;
end

if nargin < 10
    HILBERTFLAG = false;
end

if nargin < 11
    CANUSEGPU = false;
end

%%
nfft2 = nfft/2;
Hbu = Hbu(1:nfft2);
nIter = 5; %ceil(nTapers/2);
npts = pull(S,"npts");
nOverlap = round(totN*overlapPrcnt);
increment = totN-nOverlap;
nTimeWindows = 1+floor((npts(1)-totN)/increment);
lS = length(S);
AmpSmooth = NaN(nfft2,nTimeWindows,lS);
PhaseSmooth = AmpSmooth;

%%
if CANUSEGPU
    dpssSeq = gpuArray(dpssSeq);
end
i = 1;
d = pull(S(i));
[D,variance] = raw_eigen_coeffs(d,dpssSeq,nfft,nfft2,totN,...
    overlapPrcnt,detrendFlag,GAINCORRECTION,HILBERTFLAG);
[W_,Amp_,Phase_] = spec_est(D,lambda,variance,nIter);
%figure(); imagesc(squeeze(W_(:,10,:))); zoom on; colorbar;
AmpSmooth(:,:,i) = Amp_;
PhaseSmooth(:,:,i) = Phase_;

if lS < 2
    CrossSpec = [];
    return;
end

%%
nTapers = size(dpssSeq,3);
Dprev = conj(D);
Wmain = W_;
conjD = NaN(nfft2,nTimeWindows*lS,nTapers);
indivW = conjD;
for j = i+1:lS
    % get ALL autocorrelations and weights
    S_ = S(j);
    d = pull(S_);
    D = raw_eigen_coeffs(d,dpssSeq,nfft,nfft2,totN,...
        overlapPrcnt,detrendFlag,GAINCORRECTION,HILBERTFLAG);
    Amp_ = sum(Wmain.*(D.*conj(D)),3);
    Phase_ = angle(D);
    meanx = sum(Wmain.*cos(Phase_),3);
    meany = sum(Wmain.*sin(Phase_),3);
    Phase_ = atan2(meany,meanx);
    AmpSmooth(:,:,j) = Amp_;
    PhaseSmooth(:,:,j) = Phase_;
    si = (j-1)*nTimeWindows + 1;
    ei = j*nTimeWindows;
    conjD(:,si:ei,:) = conj(D);     % unweighted
    indivW(:,si:ei,:) = Wmain;
    fprintf("autocorrelation: %d\n",j);
end

if lS < 3
    CrossSpec = sum(Wmain.*(D.*Dprev),3);
    CrossSpec = conj(CrossSpec); %this conj is not a typo
    return;
end

indivW(:,1:nTimeWindows,:) = Wmain;
conjD(:,1:nTimeWindows,:) = Dprev;
refBlock = conj(conjD);         % unweighted
conjD = indivW.*conjD;          % apply weights to all channels
refBlock = circshift(refBlock,-nTimeWindows,2);
[numShifts,ncombos,excess,shuffleVec,flipVec] = shuffleVector(lS);
fprintf("getting cross spectra\n");
if numShifts < 2
    %% n == 3 (numShifts == 1)
    fprintf("num_shifts == 1\n");
    CrossSpec = conjD.*refBlock;
    CrossSpec = Hbu.*sum(CrossSpec,3);
    CrossSpec = reshape(CrossSpec,nfft2,nTimeWindows,lS);
    CrossSpec(:,:,shuffleVec) = CrossSpec;
    CrossSpec(:,:,flipVec) = conj(CrossSpec(:,:,flipVec));
    return;
end

CrossSpec = NaN(nfft2,nTimeWindows,ncombos);
ei = lS*(1:numShifts)';
si = 1 + lS*(0:numShifts-1)';
for i = 1:numShifts-1
    ei_ = ei(i);
    si_ = si(i);
    CrossSpec_ = conjD.*refBlock;
    CrossSpec_ = Hbu.*sum(CrossSpec_,3);
    CrossSpec_ = reshape(CrossSpec_,nfft2,nTimeWindows,lS);
    CrossSpec(:,:,si_:ei_) = CrossSpec_;
    refBlock = circshift(refBlock,-nTimeWindows,2);    % ready for next iteration
    fprintf("%d %d\n",i,numShifts);
end

si_ = si(numShifts);
ei_ = ei(numShifts);
if excess
    ei_ = ei_ - excess;
    conjD = conjD(:,1:nTimeWindows*(lS-excess));
    refBlock = refBlock(:,1:nTimeWindows*(lS-excess));
end
CrossSpec_ = conjD.*refBlock;
CrossSpec_ = Hbu.*sum(CrossSpec_,3);
CrossSpec_ = reshape(CrossSpec_,nfft2,nTimeWindows,[]);
CrossSpec(:,:,si_:ei_) = CrossSpec_;
CrossSpec(:,:,shuffleVec) = CrossSpec;
CrossSpec(:,:,flipVec) = conj(CrossSpec(:,:,flipVec));
end

function [D,variance] = raw_eigen_coeffs(d,dpssSeq,nfft,nfft2,totN,...
    overlapPrcnt,detrendFlag,GAINCORRECTION,HILBERTFLAG)
dcut = cutWindows(d,totN,overlapPrcnt,detrendFlag);
variance = rssq(dcut,1); %timedomain est.
dcut = dpssSeq.*dcut;
if GAINCORRECTION
    dcut = dcut./variance; %divide by root sum of squares (L2 norm)
end
variance = variance.^2;
variance = variance/(totN-1);

D = fft(dcut,nfft,1);
D = D(1:nfft2,:,:);
if ~HILBERTFLAG
    return;
end
D = -1j*D;
end

function [weights,varargout] = spec_est(D,lambda,variance,nIter)
Amp = D.*conj(D);
[weights,Amp] = get_weights(Amp,lambda,variance,nIter);
if nargout < 2
    return;
end
varargout{1} = Amp;
Phase = angle(D);
meanx = sum(weights.*cos(Phase),3);
meany = sum(weights.*sin(Phase),3);
Phase = atan2(meany,meanx);
varargout{2} = Phase;
end