function varargout = loadSubspaceCatalog(dirName,variableNames,...
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
    variableNames = [];
end

if isempty(variableNames)
    variableNames = {"tabs";"energyRatio";"z2p";"p2rms";"Neff";...
        "kurt";"skew";"z2pSynth";"synthAmpRatio";"obsAmpRatio";"ReconstructCoefficients"};
end

mainDir = fullfile("~","masa","subspace_detector",dirName);
if ~exist(mainDir,'dir')
    fprintf(2,"Directory %s doesnt exist. Are you sure you have the right name?\n",mainDir);
    return;
end

%%
cd(mainDir);
files = dir('daySubspaceSearch*mat');
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
    if ismember(thisName,["synthAmpRatio";"obsAmpRatio";"ReconstructCoefficients"])
        theseData = horzcat(S.(thisName));
        varargout{i} = theseData';
    else
        varargout{i} = vertcat(S.(thisName));
    end
end

if badN > 0
    fprintf("number of bad files: %d\n",badN);
    fNameBad = fNameBad(1:badN);
    disp(fNameBad)
end