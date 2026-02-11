function sac2txt(S,outFile)
if nargin < 2; outFile = '~/text.txt'; end

lS = length(S);
if lS > 1
    disp('input more than one station. using only the first.');
    S = S(1);
end
t = getTimeVec(S);
d = S.d;

t = datenum(t);
t = t-693960; %695422;

formatSpec = '%f %f';
str = compose(formatSpec,t,d);
str = string(str);

fileID = fopen(outFile,'w');
fprintf(fileID,'%s\n',str);
fclose(fileID);
