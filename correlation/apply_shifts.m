function shifted_data = apply_shifts(indiv_events,raw_shifts)
[winlen,nevents] = size(indiv_events);
shifted_data = zeros(winlen,nevents);

lraw = length(raw_shifts);
if lraw ~= nevents
    fprintf('dimension mismatch: columns of indiv_events not same as length of raw_shifts\n');
    return;
end

%%
for i = 1:length(raw_shifts)
    shifti = raw_shifts(i);
    if shifti > 0
        tmp = indiv_events(1:end-shifti,i);
        shifted_data(shifti+1:end,i) = tmp;
    elseif shifti < 0
        lcut = winlen + shifti;
        tmp = indiv_events(1-shifti:end,i);
        shifted_data(1:lcut,i) = tmp;
    else
        shifted_data(:,i) = indiv_events(:,i);
    end
end