function plotParticleMotion(S,cmp1,cmp2,ax)
if nargin < 4; ax = newplot; end
if nargin < 2; cmp1 = 1; end
if nargin < 3; cmp2 = 2; end

e = S(cmp1).d;
n = S(cmp2).d;
plot(ax,e,n);