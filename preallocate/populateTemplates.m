function populateTemplates(E,days,yyyy,radThresh,magThresh,base)
% Define defaults
if nargin < 6
    base = '/auto/proj/stephenh/chiles/';
end
if nargin < 5
    magThresh = 5.5;
end
if nargin < 4
    radThresh = 15;
end
if nargin < 3
    yyyy = 2014;
end
if nargin < 2
    disp('Error. Need at least two arguments in: E & days')
    return
end

nDays = length(days);
refDays = datenum(yyyy,01,00);  % A scalar
refDays = refDays + days;       % Possibly a vector

cd(base)
refla = 0.79227;
reflo = -77.9462;
for i = 1:length(E)
    origabstime(i) = datenum(E(i).origyyyy,E(i).origmm,E(i).origdd,E(i).orighh,E(i).origmmm,E(i).origsec);
    origmag(i) = E(i).origmag;
    origlat(i) = E(i).origlat;
    origlon(i) = E(i).origlon;
end

%distance from CHL1
d = distance(refla,reflo,origlat,origlon,[6378.137 0.081819190]);
dI = d <= radThresh & origmag <= magThresh;

yyyy = num2str(yyyy);
for i = 1:nDays
    doy = days(i);
    doy = num2str(doy);
    disp(['Processing: ',strcat(yyyy,doy)])
    populateDataDirectories(base,str2double(yyyy),str2double(doy));
    
    % CD into directory; determine if template structure exists
    cd(yyyy)
    cd(doy)
    masterStructure = strcat('ecua_',yyyy,'_',doy);
    
    tI = origabstime >= refDays(i) & origabstime < refDays(i)+1; % Events from just one day
    tI = dI&tI;
    if sum(tI)
        disp([num2str(sum(tI)),' events for this day.'])
        Etmp = E(tI);
        getEventTemplates(Etmp,masterStructure);
        cd ../../
    else
        disp(['No events found within specified parameters on day: ',doy])
        cd ../../
        continue
    end
end