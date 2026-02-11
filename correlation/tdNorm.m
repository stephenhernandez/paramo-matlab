function S = tdNorm(S,smoother,expo,Fs)
%
% [S,env] = tdNorm(S,smoother,expo,Fs)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Last Edited: 22 Aug. 2024

%%
if nargin < 3 
    expo = 1;
end

%% assumed the data are already filtered
isstructS = isstruct(S);
npoles = 1;
zeroPhaseFlag = true; % keep true
tw = 100; %0.008;
if isstructS
    %% get the envelope and filter
    S = detrendWaveforms(S);
    Senv = envelopeWaveforms(S);
    SenvSmooth = filterWaveforms(Senv,-inf,1/smoother,npoles,false,zeroPhaseFlag);

    %% scale each trace individually
    lS = length(S);
    for i = 1:lS
        d = S(i).d;
        d = detrend(d);
        env = SenvSmooth(i).d;

        if expo ~= 1
            env = env.^expo;
        end

        env = 1./env;
        d = taper(d.*env,tw);
        drms = rms(d);
        S(i).d = d./drms;

        [minVals,maxVals,meanVals] = minmaxmean(d);
        S(i).depmin = minVals;
        S(i).depmax = maxVals;
        S(i).depmen = meanVals;
    end
else
    if nargin < 4
        Fs = 100;
    end

    %%
    Senv = abs(hilbert(S));
    Senv = zpkFilter(Senv,-inf,1/(smoother*Fs),1,npoles,zeroPhaseFlag);
    if expo ~= 1
        Senv = Senv.^expo;
    end

    %%
    Senv = 1./Senv;
    S = taper(S.*Senv,tw);
    Srms = rms(S);
    S = S./Srms;
end