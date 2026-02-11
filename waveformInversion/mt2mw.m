function [mw,m0] = mt2mw(m)
% Get the Magnitude
mt = [m(1) m(4) m(5);
    m(4) m(2) m(6);
    m(5) m(6) m(3)];

eval3 = sort(abs(eig(mt)));


m0 = sum(eval3(2:3));
m0 = m0*1e09; %to use in conjunction with CR greens functions: 26->disp, 27->vel, 28->acc, or 21 for acc if the numbers look weird
m0 = m0/2;
mw = m02mw(m0);