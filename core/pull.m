function data = pull(S,fieldname,datatype,verboseFlag)
%
% pull pulls data like `getfield,' but more user frieldly
%
% data = pull(S,fieldname,verboseFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
if nargin < 2
    fieldname = 'd';
end

if nargin < 3
    datatype = [];
end

if nargin < 4
    verboseFlag = false;
end

%%
if verboseFlag
    disp(['Extracting: ',fieldname]);
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);

%% special cases
knownStringTypes = ["kstnm";"id";"stnm";"ntwk";"chan";"earthModel";"magType";....
    "status";"code";"knetwk";"kcmpnm";"khole";"locid";"evType";"evStatus";"evDescription";...
    "methodID";"earthModel"];
knownTimeTypes = ["ref";"t";"tStart";"tEnd";"user0"];
knownDurationTypes = ["e";"b"];

%%
if isempty(datatype)
    if ismember(fieldname,knownStringTypes)
        datatype = "";
    elseif ismember(fieldname,knownTimeTypes)
        datatype = NaT(1);
    elseif ismember(fieldname,knownDurationTypes)
        datatype = NaT(1) - NaT(1);
    else
        datatype = single(1.0);
    end
end

%%
if isstring(datatype)
    data = repmat("",lS,1);
    for i = 1:lS
        data(i) = string(S(i).(fieldname));
    end
    data = reshape(data,sizeS);
    return;
elseif isdatetime(datatype)
    var_ = S(1).(fieldname);
    if isnumeric(var_)
        data = NaN(lS,1);
    else
        data = NaT(lS,1);
    end

    data(1) = var_;
    if lS == 1
        return;
    end

    for i = 2:lS
        var_ = S(i).(fieldname);
        data(i) = var_;
    end
    data = reshape(data,sizeS);
    return;
elseif isduration(datatype)
    data = NaT(lS,1) - NaT(lS,1);
    for i = 1:lS
        data(i) = S(i).(fieldname);
    end
    data = reshape(data,sizeS);
    return; 
end

%%

if strcmp(fieldname,'d') || strcmp(fieldname,'data')
    npts = S(1).npts;
    if ~isfinite(npts)
        data = [];
        return;
    end

    % potentially recursive bit...
    npts = pull(S,'npts'); %<-- this is bad, i dont know why i did this
    npts = max(npts(:));
    data = (zeros(npts,lS)); %in one fell swoop...

    for i = 1:lS
        dtmp = S(i).(fieldname);
        nd = size(dtmp,1);

        if isfinite(S(i).npts)
            data(1:nd,i) = dtmp;
        end
    end
else
    data = (zeros(lS,1));
    for i = 1:lS
        data(i) = S(i).(fieldname);
    end
    data = reshape(data,sizeS);
end
