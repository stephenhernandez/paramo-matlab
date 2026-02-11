function [status,Sout,nSphases] = removeSphaseDuplicates(Sin)
%sort data
[Sin,chardata,datadata] = sortSin(Sin);
Sout = Sin;
nSphases = length(Sout);

if isnan(chardata{1})
    status = 0;
    nSphases = 0;
    return
else
    %find duplicates (if any)
    [~,~,IA] = unique(chardata(:,1));
    potentialRepeats = diff(IA);
    fI = find(~potentialRepeats);
    trueI = true(nSphases,1);
    
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
        %populate Sout
        chardata = chardata(trueI,:);
        datadata = datadata(trueI,:);
        lengthTrue = length(find(trueI));
        for i = 1:lengthTrue
            Sout(i).stnm = chardata{i,1};
            Sout(i).ntwk = chardata{i,2};
            Sout(i).dist = datadata(i,1);
            Sout(i).azimuth = datadata(i,2);
            Sout(i).hh = datadata(i,3);
            Sout(i).mmm = datadata(i,4);
            Sout(i).sec = datadata(i,5);
            Sout(i).res = datadata(i,6);
            Sout(i).wt = datadata(i,7);
            Sout(i).t = datadata(i,8);
        end
        Sout = Sout(1:lengthTrue);
    else
        status = 0;
    end
end