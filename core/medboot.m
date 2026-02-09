function [v,v_] = medboot(x,nboot,p)
if nargin < 3
    p = [2.5 97.5];
end

%%
[nr,nc] = size(x);

%%
if iscolumn(p)
    p = p';
end

%%
lp = length(p);
v = NaN(nc,lp);
G = randi(nr,nr,nboot);
for i = 1:nc
    x_ = x(:,i);
    x_ = x_(G);
    
    %%
    v_ = median(x_);
    v(i,:) = prctile(v_,p); % mad(v_,1);
end