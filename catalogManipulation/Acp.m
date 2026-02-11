function A = Acp(c,p,S,T)
if p == 1
    A = log(T+c) - log(S+c);
else
    A = (((T+c)^(1-p)) - ((S+c)^(1-p)))/(1-p);
end