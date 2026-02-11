function [doa,horvel] = doa_av_inversion(dT,Ginv)
%[doa,horvel,bel] = doa_av_inversion(dT,G)
slow = Ginv*dT;
doa = 90-atan2d(slow(2),slow(1))+180;
horSlow = rssq(slow(1:2));
horvel = 1./horSlow;
%bel = atan2d(slow(3),horSlow);