%Script modified from Nicholas Van Der Elst's code
%this code stops r ratio calculations once a global trigger has passed
%through a partiular spatial grid in between particular M1s and M2s.
%This prevents multiple far-field events sharing the same M1 and M2.
%this is for far-field events only.

%THE LOCAL INDEX MUST COINCIDE WITH DATA FROM LOCALCATALOG!
%function [R, Amp, M1, M2, T1, T2, gridR, Index1, Index2] = r_ratio(t,mag,gt,gmag,glat,glon,localIndex,spatialData)
function [R, Amp, M1, M2, T1, T2, gridR, Index1, Index2] = ...
    r_ratio(t,mag,gt,gmag,glat,glon,localIndex,spatialData)
% Read the spatial data for our local index
res = spatialData(1,1);
latmin = spatialData(1,2); 
latmax = spatialData(1,3);
lonmin = spatialData(1,4); 
lonmax = spatialData(1,5);

gridLat=latmin:res:latmax;
gridLon=lonmin:res:lonmax;

%Amplitude Constants
Elr0=1e-9;                          % minimum strain threshold
X0=log10(Elr0*20*3.5e9/(2*pi))+2;   % coefficients of threshold for use in script

%Spatial Constants
km2deg=90/10000;                    % coarse conversion factor
Dmin=400*km2deg;                    % distance threshold

%Convert to radians
glatrad=deg2rad(glat);
glonrad=deg2rad(glon);

%Initialize variables with lots of extra space
R=NaN(5e7,1);              	% time ratio, R
T1=R;                           % time to last quake before trigger
T2=R;                           % time to first quake after trigger
Amp=R;               		% log10 amplitude of trigger (nanometers)
Index1=R;                       % index to first quake before trigger
Index2=R;                       % index to first quake after trigger
Gindex=R;                       % index to trigger

%single=ones(length(gridLat)-1,length(gridLon)-1);
gridR=zeros(length(gridLat)-1,length(gridLon)-1);
n = 1;

%Calculate R ratio
for i=1:length(gridLat)-1
    for j=1:length(gridLon)-1
        %Get local catalog from pre-indexed matrix
        Ibin=localIndex{i,j};
        
        %Find potential triggers and calculate their Amplitudes
        D=globdist(deg2rad(gridLat(i)+res/2), deg2rad(gridLon(j)+res/2), ...
            glatrad, glonrad);
        Idist = (D >= deg2rad(Dmin));
        Ithresh = (gmag > X0+1.66*log10(rad2deg(D)));
        goodCandidates = (Idist & Ithresh);
        amp = gmag(goodCandidates)-1.66*log10(rad2deg(D(goodCandidates)))-2;
        
        %Create bin-specific global and local catalogs
        VEL = 4.3;
        %VEL = 299792.458; %<-- no longer the limit
        surfTravel = ((rad2deg(D(goodCandidates))/km2deg)/VEL)/86400;
        gtsub = gt(goodCandidates) + surfTravel;
        tsub = t(Ibin);
        goodList=find(goodCandidates);  % make an index out of logical mask
        
        %Compute R ratio
        Rbin=NaN(size(gtsub));        
        for itrig=1:length(gtsub)

	    %Select the first events before and after the arrival of the far field trigger energy
	    Index1sub=find(tsub<gtsub(itrig),1,'last');
            Index2sub=find(tsub>gtsub(itrig),1,'first');
            
            if ~isempty(Index1sub) && ~isempty(Index2sub) && Index1sub*Index2sub > 0 && Index2sub <= length(tsub) %&& Index1sub >= single(i, j) %&& T1sub <= minT && T2sub <= minT
%                single(i,j) = Index2sub;
%                Gindex(n) = goodList(itrig);
                Amp(n) = amp(itrig);
                Index1(n) = Ibin(Index1sub);
                Index2(n) = Ibin(Index2sub);
                T1(n) = gtsub(itrig)-tsub(Index1sub);%T1sub;
                T2(n) = tsub(Index2sub)-gtsub(itrig);%T2sub;
                R(n) = T2(n)./(T1(n)+T2(n));
                
                Rbin(itrig)=R(n);
                n=n+1;
            end
        end     %Processing for this spatial bin has ended.
        gridR(i,j) = mean(Rbin);
        clear Rbin
    end
end

nans=isnan(R);
R=R(~nans);
T1=T1(~nans);
T2=T2(~nans);
Index1=Index1(~nans);
Index2=Index2(~nans);
M1=mag(Index1);
M2=mag(Index2);
Amp=Amp(~nans);
