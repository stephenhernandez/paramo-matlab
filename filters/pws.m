function [phase_weighted_stack,w,linear_stack,indiv_events,phase_t] = pws(indiv_events,...
    normalizeFlag,detrendFlag,nu,Nsmooth)
if nargin < 2; normalizeFlag = true; end
if nargin < 3; detrendFlag = true; end
if nargin < 4; nu = 1; end
if nargin < 5; Nsmooth = 0; end

%%
if detrendFlag
    indiv_events = detrend(indiv_events);
end

if normalizeFlag
    indiv_events = normalizeWaveforms(indiv_events);
end

%%
if Nsmooth > 0
    box = ones(Nsmooth,1)/Nsmooth;
    linear_stack = flipud(convn(indiv_events',box));
    linear_stack = flipud(linear_stack(Nsmooth:end,:))';
else
    linear_stack = mean(indiv_events,2,"omitnan");
end

%%
if isreal(indiv_events)
    phase_t = hilbert(indiv_events);
else
    phase_t = indiv_events;
end
phase_t = angle(phase_t);           %phase in radians
phase_t = exp(1j*phase_t);

%%
if Nsmooth > 0
    phase_t = flipud(convn(phase_t',box));
    w = abs(flipud(phase_t(Nsmooth:end,:))');
else
    w = abs(mean(phase_t,2,"omitnan")).^nu; %mean phase
end

%%
phase_weighted_stack = linear_stack.*w;