function [S,Hd] = filterWaveforms(S,varargin)
%
% filterWaveforms return structure with filtered waveforms
%
% S = filterWaveforms(S,lfc,hfc,npoles,tw,zeroPhaseFlag,diffFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 1
    fprintf(2,'not enough inputs, returning empty structure\n')
    S = populateWaveforms();
    return;
end

%%
nVarargin = length(varargin);
functionDefaults = {...
    -inf,...    % lfc
    -inf,...    % hfc
    4,...       % npoles
    false,...   % taperWidth
    false,...   % zeroPhaseFlag
    false};     % differentiateFlag

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[lfc,hfc,npoles,tw,zeroPhaseFlag,diffFlag] = deal(optsToUse{:});

%%
if isempty(lfc)
    lfc = -inf;
end

if isempty(hfc)
    hfc = -inf;
end

if isempty(tw)
    tw = false;
end

if isempty(npoles)
    npoles = 4;
end

%%
sizeS = size(S);
S = S(:);

if diffFlag
    S = differentiateWaveforms(S);
end

%%
cornersfin = isfinite([lfc hfc]);
if ~any(cornersfin)
    fprintf(2,"no valid corners input, doing nothing\n");
    return;
end

%%
deltas = pull(S,"delta");
Fs = 1./deltas;
fI = Fs >= 1;
Fs(fI) = round(Fs(fI));
uniqFs = unique(Fs);
luniqFs = sum(isfinite(uniqFs));

%%
for i = 1:luniqFs
    uniqFs_ = uniqFs(i);
    Hd = [];
    if ~zeroPhaseFlag
        Hd = zpkOperator(lfc,hfc,uniqFs_,npoles);
    else
        [~,sos,g] = zpkOperator(lfc,hfc,uniqFs_,npoles);
    end

    %%
    fI = Fs == uniqFs_;
    S_ = S(fI);                             %waveforms with a uniq Fs
    npts = pull(S_,"npts");
    uniqnpts = unique(npts);
    luniqnpts = length(uniqnpts);
    for j = 1:luniqnpts
        uniqnpts_ = uniqnpts(j);
        nI = uniqnpts_ == npts;
        nn = sum(nI);
        S__ = S_(nI);

        %% high-pass filtering, remove trend
        if cornersfin(1)
            S__ = demeanWaveforms(S__); %ignores NaNs
        end

        %% re-populate with filtered versions
        d_filtered = double(pull(S__));
        for k = 1:nn
            d_tmp = d_filtered(:,k);
            goodI = isfinite(d_tmp);
            d_tmp = d_tmp(goodI);

            %% filter here (with optional zero-phase filtering (two pass))
            if ~zeroPhaseFlag
                d_tmp = filter(Hd,d_tmp);
            else
                warning off signal:filtfilt:ParseSOS
                try
                    d_tmp = filtfilt(sos,g,d_tmp);
                catch
                    fprintf("not enough data to filter, trying with zeros...\n");
                    d_tmp = d_filtered(:,k);
                    badI = ~isfinite(d_tmp);
                    d_tmp(badI) = 0;
                    d_tmp = filtfilt(sos,g,d_tmp);
                    d_tmp = d_tmp(goodI);
                end
            end
            d_filtered(goodI,k) = d_tmp; %fill with actual filtered data

            %%
            S__(k).d = d_filtered(:,k);
            [minVals,maxVals,meanVals] = minmaxmean(d_filtered(:,k));
            S__(k).depmin = minVals;
            S__(k).depmax = maxVals;
            S__(k).depmen = meanVals;
        end
        S_(nI) = S__;
    end
    S(fI) = S_;
end

%% taper if desired
if tw
    S = taperWaveforms(S,tw);
end
S = reshape(S,sizeS);