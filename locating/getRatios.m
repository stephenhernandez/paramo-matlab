function obsRatio = getRatios(amps,divFlag)
if nargin < 2
    divFlag = ~true;
end

%%
n = 0;
[Namp,dim2] = size(amps);
Np = 0.5*Namp*(Namp-1);
ampTmp = amps(2:end,1:dim2);
Nttmp = size(ampTmp,1);
obsRatio= zeros(Np,dim2);

%%
if divFlag
    for i = 1:Namp-1
        for j = 1:Nttmp
            n = n+1;
            obsRatio(n,1:dim2) = amps(i,1:dim2)./ampTmp(j,1:dim2);
        end
        ampTmp = ampTmp(2:end,1:dim2);
        Nttmp = size(ampTmp,1);
    end
else
    for i = 1:Namp-1
        for j = 1:Nttmp
            n = n+1;
            obsRatio(n,1:dim2) = amps(i,1:dim2)-ampTmp(j,1:dim2);
        end
        ampTmp = ampTmp(2:end,1:dim2);
        Nttmp = size(ampTmp,1);
    end
end