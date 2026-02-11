function generateHtmlTabs(tStart,tEnd)
if nargin < 2; tStart = dn2dt(now); end
if nargin < 1; tStart = dn2dt(now); end

regions = "cotopaxi";
lRegions = length(regions);
days = tStart:tEnd;
lDays = length(days);

%% load seiscomp catalog
[tOrig,eqlatOrig,eqlonOrig,eqdepthOrig,eqmagOrig,idOrig,rmsOrig,azgapOrig,~,nMLvOrig,...
    ~,~,~,~,~,~,~,~,~,eqTypeOrig,creationTimeOrig] = readCat1();

%% define filters
minMag = 1;
maxRMS = 2;
maxAzGap = 180;
tWinStart = [-365 -28 -7 0];
tWinEnd = [0 0 0 1];

%% loop
for i = 1:lRegions
    regionName = regions(i);
    boundaryBox = getRegionSpatialDimensions(regionName);
    minLon = boundaryBox(1);
    maxLon = boundaryBox(2);
    minLat = boundaryBox(3);
    maxLat = boundaryBox(4);
    
    for j = 1:length(tWinStart)
        data = [];
        tWinStart_ = tWinStart(j);
        tWinEnd_ = tWinEnd(j);
        for k = 1:lDays
            currentDay = days(k);
            [yyyy,mmm,ddd] = datevec(currentDay);
            dirName = char(strcat("~/public_html/",regions(i),"/",num2str(yyyy),"/",num2str(mmm),"/",num2str(ddd),"/"));
            if ~exist(dirName,'dir')
                mkdir(dirName);
            end
            
            %% filter catalog
            tI = eqlatOrig >= minLat & eqlatOrig <= maxLat & eqlonOrig >= minLon ...
                & eqlonOrig <= maxLon & rmsOrig <= maxRMS & eqmagOrig >= minMag ...
                & (azgapOrig <= maxAzGap | nMLvOrig >= minNMLV) & ...
                tOrig >= currentDay+tWinStart_ & tOrig >= currentDay+tWinEnd_;
            if sum(tI)
                t = tOrig(tI);
                eqlat = eqlatOrig(tI);
                eqlon = eqlonOrig(tI);
                eqdepth = eqdepthOrig(tI);
                eqmag = eqmagOrig(tI);
                %ids = idOrig(tI);
                %rmssec = rmsOrig(tI);
                %azgap = azgapOrig(tI);
                %eqType = eqTypeOrig(tI);
                %creation = creationTimeOrig(tI);
                data_ = [t eqlat eqlon eqdepth eqmag];
                data = cat(3,data,data_);
            end
        end
        
        %% plot data
        % ive generated the data, now need to pass it to following function
        % so that it can load figure (once) and loop through that data 3D
        % matrix and progressively plot
        plotEcuadorSeismicity(tRef,data,boundaryBox,urbanFlag,printFlag)
    end
end