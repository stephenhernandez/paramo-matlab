function [S,w,linear_stack] = pwsWaveforms(S,normalizeFlag,detrendFlag)
if nargin < 2; normalizeFlag = true; end
if nargin < 3; detrendFlag = true; end

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
isstructS = isstruct(S);
if isstructS
    npts = pull(S,'npts');
    if numel(unique(npts)) == 1
        data = pull(S);
        [stack,w,linear_stack] = pws(data,normalizeFlag,detrendFlag);
        S = S(1);
        S.d = stack;
    else
        fprintf(2,'some of the traces arent the same size, quitting\n');
        return;
    end
else
    % input data are not waveform structs. treat as matrix where columns
    % are individual traces.
    [S,w,linear_stack] = pws(S,normalizeFlag,detrendFlag);
end
