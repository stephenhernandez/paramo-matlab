function obsDiff = getDD(tt,minusOperator)
if nargin < 2
    minusOperator = true;
end

%%
n = 0;
[Ntt,dim2] = size(tt);
Np = 0.5*Ntt*(Ntt-1);
ttmp = tt(2:end,1:dim2);
Nttmp = size(ttmp,1);
obsDiff= zeros(Np,dim2);


for i = 1:Ntt-1
    for j = 1:Nttmp
        %disp([i j])
        n = n+1;
        if minusOperator
            obsDiff(n,1:dim2) = tt(i,1:dim2)-ttmp(j,1:dim2);
        else
            obsDiff(n,1:dim2) = tt(i,1:dim2)./ttmp(j,1:dim2);
        end
    end
    ttmp = ttmp(2:end,1:dim2);
    Nttmp = size(ttmp,1);
end