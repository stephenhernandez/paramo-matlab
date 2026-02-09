function S = push(S,fieldname,data,verboseFlag)
%
% push: push data to structure
%
% S = push(S,fieldname,data,verboseFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
if nargin < 2
    fieldname = "d";
end

if nargin < 4
    verboseFlag = false;
end

%%
if verboseFlag
    fprintf("Pushing: %s\n",fieldname);
end

%%
[lr,lc] = size(data);
S = S(:);
lS = length(S);
if strcmp(fieldname,"d") || strcmp(fieldname,"data")
    if lc ~= lS
        fprintf(2,"number of columns dont match number of S entries\n");
        return;
    end

    for i = 1:lS
        S(i).(fieldname) = data(:,i);
    end
else
    if lr ~= lS
        fprintf(2,"number of rows dont match number of S entries\n");
        return;
    end

    for i = 1:lS
        S(i).(fieldname) = data(i);
    end
end