function d = resampleSH(d,delta,newFs,verboseFlag,detrendFlag)
if nargin < 4
    verboseFlag = true;
end

if nargin < 5
    detrendFlag = false;
end

Fs = 1/delta;

%%
ratio = Fs/newFs;
if ratio == 1
    if verboseFlag
        fprintf('already at new Fs, no need to resample, doing nothing\n');
    end
    return;
end

%% remove trend
if detrendFlag
    t = (1:length(d))';
    p_ = polyfit(t,d,1); %fit linear trend
    dSynth = p_(1)*t + p_(2);
    d = d - dSynth;
else
    %centroidM = mean(d,"omitnan");
    centroidM = median(d,"omitnan");
    d = d - centroidM;
end

%%
modulus = mod(Fs,newFs);
if ~modulus
    % clean downsampling
    ratio = round(ratio);
    facts = factor(ratio);
    facts = fliplr(facts);

    for j = 1:length(facts)
        downFact = facts(j);
        d = resample(d,1,downFact);
    end

    %%
    if verboseFlag
        fprintf('decimated data to %u\n',newFs);
    end
else
    % either upsample or non-clean downsample
    origP = newFs/Fs;
    [newP,Q] = rat(origP,1e-5);

    try
        %floor(2^31/Q);
        d = resample(d,newP,Q);
    catch
        newP = floor(2^31/Q);
        d = resample(d,newP,Q);
        % if ~verboseFlag
        %     fprintf('gyatt damn.. something done gone wrong\n');
        % end
        % disp(['ID: ' ME.identifier])
        % rethrow(ME);
    end
end

%%
if detrendFlag
    t = (1:length(d))';
    dSynth = p_(1)*t + p_(2);
    d = d + dSynth; %add back the trend you previously removed
else
    d = d + centroidM;
end
