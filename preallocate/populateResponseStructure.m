function responseStructure = populateResponseStructure(sizeR)

responseStructure.Tstart = NaT;
responseStructure.Tend = NaT;
responseStructure.Stel = NaN;
responseStructure.Stlo = NaN;
responseStructure.Stla = NaN;

responseStructure.knetwk = "";
responseStructure.kstnm = "";
responseStructure.khole = "";
responseStructure.kcmpnm = "";
responseStructure.instType = "";

responseStructure.N = 0;
responseStructure.H = [];

responseStructure.zeros = [];
responseStructure.poles = [];
responseStructure.constant = NaN;
responseStructure.gain =NaN;
responseStructure.sensitivity = NaN;
responseStructure.A0 = NaN;
responseStructure.Fs = NaN;

%%
if numel(sizeR) == 1
    responseStructure = repmat(responseStructure,sizeR,1); %make column vector
else
    responseStructure = repmat(responseStructure,sizeR);
end