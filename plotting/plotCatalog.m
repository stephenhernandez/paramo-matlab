function ax = plotEcuadorianCatalog(varargin)
% Example: [] =
% generateSeismicityAnimationFrames(catalogData,boundaryBox,tStart,tEnd,regionName,demName,minMag,
% depthPanel,cumPanel,maxDepth,addSimbologia,allFrames,hillShadeFlag,magFact,fontSize,fstring,minLevInc);

% define defaults
functionDefaults = {...
    catalogData,...           	% catalogData = [t,eqlat,eqlon,eqdepth,eqmag,id]
    boundaryBox,...             % boundaryBox
    dn2dt(floor(now)),...       % start time
    dn2dt(ceil(now)),...        % end time
    'chiles',...                % regionName
    [],...                      % boundaryBox
    true,...                    % SC3Flag
    2,...                       % min mag
    false,...                   % plot depth panel?
    true,...                    % plot cum panel?
    100,...                     % max depth
    true,...                    % add legend?
    true,...                    % all frames?
    true,...                    % hillshading?
    20,...                      % mag amplification factor
    18,...                      % text font size
    '-djpeg',...                % output format
    [],...                      % minLevInc
    false};                     % diasFlag

% deal variables
optsToUse = functionDefaults;
nVarargin = length(varargin);
optsToUse(1:nVarargin) = varargin;
[catalogData,boundaryBox,eqmag,tStart,tEnd,regionName,boundaryBox,SC3Flag,...
    minMag,depthPanel,cumPanel,maxDepth,...
    addSimbologia,allFrames,hillShadeFlag,magFact,fontSize,fstring,...
    minLevInc,diasFlag] = deal(optsToUse{:});