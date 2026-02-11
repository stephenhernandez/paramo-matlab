function G = Gvdcc(n)
%make the vandecar + crosson sparse matrix
N = n-1;
N_ = N:-1:1;
rows = n*N/2;
r1 = rows+1;
mR = 1:rows;
fC = mR;
fC2 = fC;
lN = length(N_);
cs = 2:n;
rs = 1+cumsum(N_);
rs = circshift(rs,[0 1]);
rs(1) = 1;

for i = 1:lN
    fC(rs(i):rs(i)+N_(i)-1) = cs(i):cs(i)+N_(i)-1;
    fC2(rs(i):rs(i)+N_(i)-1) = i*ones(1,N_(i));
end
R2 = r1*ones(n,1);
Clast = 1:n;
G1 = sparse(mR,fC,-ones(rows,1),r1,n);
G2 = sparse(mR,fC2,ones(rows,1),r1,n);
G3 = sparse(R2,Clast,ones(n,1),r1,n);
G = G1+G2+G3;