function tvec = getTimeVec(S,dnFlag)
%
% getTimeVec returns vector of datetime objects with start, end, and
% sampling rate as specified in input structure
%
% tvec = getTimeVec(S,dnFlag)
% dnFlag: optional (default = false) conversion to datenum format (floats)
%
% only one time vector returned, that associated with first input structure
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    dnFlag = false;
end

%%
S = S(1);
ref = S.ref;
if isnat(ref)
    disp('ref is not a datetime object.');
    tvec = [];
    return;
end

tvec = (S.ref:seconds(S.delta):S.ref+S.e)';
if dnFlag
    tvec = datenum(tvec);
end
tvec = tvec(1:S.npts);
