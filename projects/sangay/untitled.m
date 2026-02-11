% sangay VASR
clear; 

cd ~/research/now/sangay/

%%
load SAGA_VASR.mat
z2p2(t2 <= datetime(2018,11,11)) = z2p2(t2 <= datetime(2018,11,11))*4;
close all; zI = z2p2 >= 1e4; figure(); plot(t2(zI),1:sum(zI),'.'); zoom on;
t3 = t2(zI);
z2p3 = z2p2(zI);
tI = ismembertol(seconds(t3-min([min(t3) min(tBDF)])),seconds(tBDF-min([min(t3) min(tBDF)])),1,'DataScale',1);
sum(tI)
t4 = t3(tI);
z2p4 = z2p3(tI);
tI = ismembertol(seconds(tBDF-min([min(t4) min(tBDF)])),seconds(t4-min([min(t4) min(tBDF)])),1,'DataScale',1);
sum(tI)
z2pBDF2 = z2pBDF(tI);
VASR = z2pBDF2./z2p4;
figure(); semilogy(t4,VASR,'.'); zoom on; grid on;