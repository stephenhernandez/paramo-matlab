%april23_2019
clear; close all; clc;
yyyy = 2018;
iStart = 160; 
iEnd = 169; 
newFs = 16; 
ampFact = 50; 
cc = []; 
caus = cc; 
acaus = cc;
net = "9D";
for i = iStart:iEnd %this one was with weights based on hf
    
    S(1) = loadSacData(datetime(yyyy,01,i),1,"SN11","HHZ",net);
    S(2,1) = loadSacData(datetime(yyyy,01,i),1,"SN12","HHZ",net);
    [San(i-iStart+1),tdum,cc_,caus_,acaus_] = stackNoiseCorr(S,newFs,1,2,2^16,0.25,1,0,false,false);
    cc = [cc cc_]; 
    caus = [caus caus_]; 
    acaus = [acaus acaus_];
end

data = pull(San);
tmp_stack(:,1) = plot_family(data,1:size(data,2),ampFact,newFs);
p2rms = peak2rms(cc)'; pI = p2rms >= 5;

data = cc(:,pI); 
tmp_stack(:,2) = plot_family(data,1:size(data,2),ampFact,newFs);

figure(); imagesc(cc(:,pI)'); axis xy; colorbar; zoom on;
figure(); imagesc(normalize_traces(cumsum(cc(:,pI)')')'); axis xy; colorbar; zoom on;

data = normalize_traces(cumsum(cc(:,pI)')'); 
tmp_stack(:,3) = plot_family(data,1:size(data,2),ampFact,newFs);

w1 = p2rms / sum(p2rms);
clear stack
stack(:,1) = mean(cc,2);
stack(:,2) = sum(w1.*(cc'))';
figure(); plot(normalize_traces(stack),'.'); zoom on;
sum(pI)
peak2rms(stack)