function [c,p,n,P,cvec,pvec,logger] = molmle(t,tStart,tStop,C,P)
% p: omori decay exponent
% c: value in seconds since mainshock
% t, tStart, tStop: units should be in seconds since mainshock

tI = t>= tStart & t <= tStop;
n = sum(tI);
t = t(tI);
Nlog = 351;
pvec = linspace(-0.2,P,Nlog)';
cvec = logspace(0,log10(C),Nlog)';

logger = NaN.*cvec';
P = 0.*meshgrid(cvec,pvec);

for j = 1:length(cvec)
    c_ = cvec(j);
    logger(j) = sum(log(1+t./c_));
end

%figure(); plot(logger,'.'); zoom on; 
for i = 1:length(pvec)
    p_ = pvec(i);
    D_ = Dcp(cvec,p_,tStart,tStop)';
    P(i,:) = n*log(D_) - p_ .* logger;
end

[~,maxCol] = max(max(P));
[~,maxRow] = max(max(P,[],2));

c = cvec(maxCol);
p = pvec(maxRow);

P = P/n;
P = P - min(min(P));
P = exp(P);
P = P/sum(sum(P));

%%
% % % % function [c,p,n,P,cvecOrig,pvecOrig] = molmle(t2,S,T,maxC,maxP)
% % % % % p: omori decay exponent
% % % % % c: value in seconds since mainshock
% % % % % t, tStart, tStop: units should be in seconds since mainshock
% % % % 
% % % % tI = t2 >= S & t2 <= T;
% % % % n = sum(tI);
% % % % t2 = t2(tI);
% % % % NP = 101;
% % % % NC = 201;
% % % % pvecOrig = linspace(0.2,maxP,10*NP)';
% % % % cvecOrig = linspace(1,(maxC),10*NC)';
% % % % 
% % % % logger = NaN.*cvecOrig';            % preallocate
% % % % [cvec,pvec] = meshgrid(cvecOrig,pvecOrig);   % preallocate
% % % % origShape = size(cvec);
% % % % cvec = cvec(:);
% % % % pvec = pvec(:);
% % % % 
% % % % %x-coordinate/left-right/nCols == c
% % % % %y-coordinate/up-down/nRows == p
% % % % P = NaN*cvec;
% % % % lP = length(P);
% % % % for i = 1:lP
% % % %     c_ = cvec(i);
% % % %     p_ = pvec(i);
% % % %     logger = sum(log(1+t2./c_));
% % % %     
% % % %     D_ = Dcp(c_,p_,S,T)';
% % % %     P(i) = -n*log(D_) + p_ .* logger;
% % % % end
% % % % 
% % % % Porig = P;
% % % % [minValue,mI] = min(P,[],1,"omitnan")
% % % % maxValue = max(P,[],1,"omitnan")
% % % % c = cvec(mI)
% % % % p = pvec(mI)
% % % % 
% % % % %%
% % % % %P = P/n;
% % % % %P = P - minValue;
% % % % P = exp(P);
% % % % P = P/sum(sum(P));
% % % % P = reshape(P,origShape);