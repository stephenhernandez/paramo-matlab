function v = stdboot(x,nboot)
[nr,nc] = size(x);

%%
v = NaN(nc,1);
G = randi(nr,nr,nboot);
for i = 1:nc
    x_ = x(:,i);
    x_ = x_(G);
    
    %%
    v_ = var(x_);
    v(i) = sqrt(mean(v_));
end