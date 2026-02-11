function [dc,clvd] = percentDC(m)

eigenvalues = eig([m(1) m(4) m(5); 
                   m(4) m(2) m(6); 
                   m(5) m(6) m(3)]);
mineig = min(abs(eigenvalues));
maxeig = max(abs(eigenvalues));
clvd = 200*mineig/maxeig;
dc = 100-clvd;