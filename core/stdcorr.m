function [v,G] = stdcorr(x)
[nr,nc] = size(x);

%%
v = NaN(nc,1);
G = zeros(nr);
for i = 1:nc
    x_ = x(:,i);
    c = xcov(x_);
    c = c(nr:end);
    
    %%
    G(:,1) = c;
    
    %%
    for j = 2:nr
        G(j:end,j) = c(1:end-j+1);
        G(j,j:end) = c(1:end-j+1);
    end
    
    %%
    v(i) = sqrt(mean(mean(G)));
end