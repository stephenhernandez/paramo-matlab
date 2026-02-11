function [tFake,cumL,tNew,iFake,lambda] = molrnd(N,c,p,tStart,tStop)
% c,tStart,tStop: in seconds

dx = 1;
D = Dcp(c(end),p(end),tStart,tStop);
tNew = (tStart:dx:tStop)';
lambda = D./((1+(tNew./c(end))).^p(end));
cumL = cumtrapz(lambda);

iFake = sort(rand(N,1));
tFake = interp1(cumL,tNew,iFake);