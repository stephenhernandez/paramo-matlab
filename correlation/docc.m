function [maxccp,plags,maxccn,nlags] = docc(data,verboseFlag)

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2; verboseFlag = false; end

[winlen,n] = size(data);
ncol = n*(n-1)*0.5;
disp(['Window Length: ',num2str(winlen)]);
disp(['Number of Events: ',num2str(n)]);
disp(['Number of Combos: ',num2str(ncol)]);

%%
ei = fliplr(1:n-1);
ei = cumsum(ei);
si = ei+1;
si = [1 si];
si(end) = [];

%%
maxccp = NaN(ncol,1);
plags = maxccp;
nlags = maxccp;
maxccn = maxccp;

%%
data = normalizeWaveforms(data);
refBlock = flipud(data);

%% get the job done
if verboseFlag
    for i = 1:n-1
        disp(i);
        si_ = si(i);
        ei_ = ei(i);
        master = data(:,i);
        refBlock = refBlock(:,2:end);
        Ctmp = convn(master,refBlock);
        
        % positive lags
        [maxccp_,plags_] = max(Ctmp);
        maxccp(si_:ei_) = maxccp_;
        plags(si_:ei_) = plags_;
        
        % negative lags
        [maxccn_,nlags_] = min(Ctmp);
        maxccn(si_:ei_) = -maxccn_;
        nlags(si_:ei_) = nlags_;
    end
else
    for i = 1:n-1
        si_ = si(i);
        ei_ = ei(i);
        master = data(:,i);
        refBlock = refBlock(:,2:end);
        Ctmp = convn(master,refBlock);
        
        % positive lags
        [maxccp_,plags_] = max(Ctmp);
        maxccp(si_:ei_) = maxccp_;
        plags(si_:ei_) = plags_;
        
        % negative lags
        [maxccn_,nlags_] = min(Ctmp);
        maxccn(si_:ei_) = -maxccn_;
        nlags(si_:ei_) = nlags_;
    end
end

%%
plags = plags - winlen;
nlags = nlags - winlen;
