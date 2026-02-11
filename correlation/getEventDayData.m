function [S,newPhaseStruct] = getEventDayData(Pphases,yyyy,mm,dd,maxNStations,resThresh,maxDist,threeFlag,rawDataDir)
if nargin < 9; rawDataDir = '~/rawdata'; end
if nargin < 8; threeFlag = false; end
if nargin < 7; maxDist = 100; end

lP = length(Pphases);
i = 0;
n = 0;

%% preallocate memory
newPhaseStruct = Pphases;
S = populateWaveforms(lP);

%%
while n < maxNStations && i < lP
    i = i+1;
    kstnm_ = Pphases(i).stnm; %kstnm_ = char(kstnm_);
    ntwk_ = Pphases(i).ntwk; %ntwk_ = char(ntwk_);
    chan_ = Pphases(i).chan; %chan_ = char(chan_);
    locID_ = Pphases(i).locid; %chan_ = char(chan_);
    dist_ = 111.19*Pphases(i).dist;
    %chan_ = [chan_,'Z'];
    res_ = abs(Pphases(i).res);
    
    fprintf('searching for %s\n',kstnm_);
    
    if res_ < resThresh && dist_ <= maxDist
%         if strcmp(kstnm_,'ISPT')
%             locID_ = "00";
%         end
%         if strcmp(kstnm_,'CABP')
%             locID_ = "00";
%         end
%         if strcmp(kstnm_,'FLFR')
%             locID_ = "00";
%         end
%         if strcmp(kstnm_,'PAYG')
%             locID_ = "00";
%         end
        if threeFlag
            chanBase = char(chan_);
            chanBase = chanBase(:,1:2);
            S_ = loadWaveforms(datetime(yyyy,mm,dd),1,kstnm_,...
                [strcat(chanBase,"E"),strcat(chanBase,"N"),strcat(chanBase,"Z")],...
                [ntwk_;"EC"],locID_,true,true,rawDataDir);
        else
            S_ = loadWaveforms(datetime(yyyy,mm,dd),1,kstnm_,chan_,ntwk_,...
                locID_,0,0,rawDataDir);
        end
        %locID_ = "";
        if isnat(S_.ref)
            disp([chan_,' does not exist on this day']);
        else
            n = n+1;
            S(n) = S_;
            newPhaseStruct(n) = Pphases(i);
        end
    end
end
S = S(1:n);
newPhaseStruct = newPhaseStruct(1:n);
