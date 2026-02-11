function varargout = loadRepeaterCatalog(dirName,variableNames,...
    dayStart,dayEnd)

%
% loadRepeaterCatalog(dirName,varargin)
%

%%
if nargin < 4
    dayEnd = datetime(2050,01,01);
end

if nargin < 3
    dayStart = datetime(1990,01,01);
end

if nargin < 2
    variableNames = {"tMain";"ccMain";"ampMain";"dMag";"magMain";...
        "templateNumber";"madMain";"nUsedMain"};
end

mainDir = fullfile("~","masa","template_search",dirName);
if ~exist(mainDir,'dir')
    fprintf(2,"Directory %s doesnt exist. Are you sure you have the right name?\n",mainDir);
    return;
end

%%
cd(mainDir);
files = dir('dayTemplateSearch*mat');
lfiles = length(files);
yyyy = NaN(lfiles,1);
mm = yyyy;
dd = yyyy;
for i = 1:lfiles
    fName = files(i).name;
    fName_split = strsplit(fName,"_");
    date_ = fName_split{2};
    date_split = strsplit(date_,".");
    yyyy(i) = str2double(date_split{1});
    mm(i) = str2double(date_split{2});
    dd(i) = str2double(date_split{3});
end

tFile = datetime(yyyy,mm,dd);
goodI = tFile >= dayStart & tFile <= dayEnd;
if ~sum(goodI)
    fprintf("no files match desired time ragne of: %s - %s\n",dayStart,dayEnd);
    return;
end

files(~goodI) = [];
lfiles = length(files);
varStrings =  string(variableNames);
lNames = length(varStrings);
if nargout ~= lNames
    fprintf(2,"mismatch between number of variables in to number of variables out\n");
    return;
end

% preallocate
S = cell(lfiles,1);
n = 0;
fNameBad = repmat("",10,1);
badN = 0;
for i = 1:lfiles
    fName = files(i).name;
    fprintf("%s\n",fName);
    try
        n = n + 1;
        S{n} = load(fName,variableNames{:});
    catch
        fprintf("Could not load file: %s\n",fName);
        n = n - 1;
        badN = badN + 1;
        fNameBad(badN) = fName;
        continue;
    end

end
S = S(1:n);
S = [S{:}];
S = S';

varargout = cell(1,nargout);
for i = 1:lNames
    thisName = varStrings(i);
    varargout{i} = cat(1,S(:).(thisName));
end

if badN > 0
    fprintf("number of bad files: %d\n",badN);
    fNameBad = fNameBad(1:badN);
    disp(fNameBad)
end

%%
% nFloatVariables = lNames;
% nMax = 1e6;
% [existT,tMainI] = ismember("tMain",varStrings);
% if existT
%     tMain = NaT(nMax,1);
%     nFloatVariables = nFloatVariables - 1;
% end
% 
% [existEvid,evidMainI] = ismember("evidMain",varStrings);
% if existEvid
%     evidMain = repmat("",nMax,1);
%     nFloatVariables = nFloatVariables - 1;
% end
% 
% if nFloatVariables
%     floaters = NaN(nMax,nFloatVariables);
% end
% 
% %get the job done
% %as is, is very slow. rewrite to separate float matrix from nat matrix from
% %string matrix
% n = 1;
% for i = 1:lfiles
%     disp(i);
%     fName = files(i).name;
%     S_ = load(fName,variableNames{:});
%     try
%         thisVarName = varStrings(1);
%         dummy = S_.(thisVarName);
%         thisN = length(dummy);
%     catch
%         fprintf(2,'variable \"%s\" does not exist in file: %s\n',varStrings(1),fName);
%         return;
%     end
% 
%     %%
%     floaterIndex = 0; %reset to 0 for each file...
%     if strcmp(thisVarName,"tMain")
%         tMain(n:n+thisN-1,1) = dummy;
%     elseif strcmp(thisVarName,"evidMain")
%         evidMain(n:n+thisN-1,1) = dummy;
%     else
%         floaterIndex = floaterIndex+1;
%         floaters(n:n+thisN-1,floaterIndex) = dummy;
%     end
%     n = n+thisN;
% 
%     %%
%     if lNames < 2
%         continue;
%     end
% 
%     %%
%     n = n-thisN; %reset n for variables to follow...
%     for j = 2:lNames
%         thisVarName = varStrings(j);
%         dummy = S_.(thisVarName);
%         if strcmp(thisVarName,"tMain")
%             tMain(n:n+thisN-1,1) = dummy;
%         elseif strcmp(thisVarName,"evidMain")
%             evidMain(n:n+thisN-1,1) = dummy;
%         else
%             floaterIndex = floaterIndex+1;
%             floaters(n:n+thisN-1,floaterIndex) = dummy;
%         end
%     end
%     n = n+thisN;
% end
% n = n-1;
% if existT
%     varargout{tMainI} = tMain(1:n);
% end
% 
% if existEvid
%     varargout{evidMainI} = evidMain(1:n);
% end
% 
% if nFloatVariables
%     floatVarNames = varStrings;
%     [lia,locb] = ismember(["tMain";"evidMain"],floatVarNames);
%     floatVarNames(locb(lia)) = [];
%     [lia,locb] = ismember(floatVarNames,varStrings);
%     if sum(lia) ~= nFloatVariables
%         fprintf(2,"unexpected number of float variables\n");
%         return;
%     end
%     locb = locb(lia);
%     for i = 1:nFloatVariables
%         varargout{locb(i)} = floaters(1:n,i);
%     end
% end