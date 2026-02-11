function [vr,vr_l,vr_u] = varRed(obs,synth,nboot)
% vr = 2*dot(obs,synth)-1;

residual = obs - synth; %synth - obs;
max_i = length(residual);
bootI = randi(max_i,max_i,nboot);
err_boot = residual(bootI);
obs_boot = obs(bootI);
vr = 100*(1-(sum(residual.^2)/sum(obs.^2)));
%vr = dot(obs,synth)/(norm(synth)*norm(obs));
vr_boot = (1-(sum(err_boot.^2)./sum(obs_boot.^2)));
vr_l = prctile(vr_boot,2.5);
vr_u = prctile(vr_boot,97.5);