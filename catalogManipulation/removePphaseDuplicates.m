function [status,Pout] = removePphaseDuplicates(Pin)

%%
% this function is completely fucked up. i dont know what happened. i
% suppose i had run into files where a stnm/ntwk combo was repeated, maybe
% the analyst picked the same SNCL twice, but ended up saving multiple
% values for the p-wave arrival, and so, i needed to pick one, so, a long
% time ago i decided to choose a pick based on its residual value (the
% smaller of all the options) and delete the rest. but it seems like a
% really clunky implementation so im going to remove it, and, when i come
% across this error in the future, i will improve the code but until then
% its not necessary to use this function
% stephen hernandez, sat. 03 aug., 2019

%%
%function [status,Pout] = removePphaseDuplicates(Pin,sortFlag)
%sort data
% if nargin < 2
%     sortFlag = false;
% end
%[Pin,chardata,datadata] = sortPin(Pin,sortFlag);
[Pin,strdata] = sortPin(Pin); %sorts according to stnm/ntwk combos
Pout = Pin;
nPphases = length(Pout);

%find duplicates (if any)
[~,~,IA] = unique(strdata); %chardata(:,1));
potentialRepeats = diff(IA);
fI = find(~potentialRepeats);
trueI = true(nPphases,1);

%remove duplicates (if any)
if any(fI)
    status = 1;
    for i = 1:length(fI)
        trueI(fI(i)) = false;
        trueI(fI(i)+1) = false;
        restmp1 = datadata(fI(i),6);
        restmp2 = datadata(fI(i)+1,6);
        if restmp1 < restmp2
            trueI(fI(i)) = true;
        else
            trueI(fI(i)+1) = true;
        end
    end
    
    %% populate Pout
    chardata = chardata(trueI,:);
    datadata = datadata(trueI,:);
    lengthTrue = length(find(trueI));
    for i = 1:lengthTrue
        Pout(i).stnm = chardata{i,1};
        Pout(i).ntwk = chardata{i,2};
        Pout(i).chan = chardata{i,3};
        Pout(i).dist = datadata(i,1);
        Pout(i).azimuth = datadata(i,2);
        Pout(i).res = datadata(i,6);
        Pout(i).wt = datadata(i,7);
        Pout(i).t = datadata(i,8);
        Pout(i).polarity = datadata(i,8);
    end
    Pout = Pout(1:lengthTrue);
else
    status = 0;
end