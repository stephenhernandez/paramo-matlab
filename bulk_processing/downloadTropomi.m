%download tropomi L3 products
function downloadTropomi(tStart)
dataHome = fullfile("~","masa","tropomi"); %tharp
cd(dataHome);

nowTime = dateshift(datetime('now') + hours(5),'start','day');
if nargin < 1
    tStart = nowTime-8;
end
tEnd = nowTime-1;
days = (tStart:tEnd)';
ldays = length(days);

Qualities = ["NRTI";"RPRO";"OFFL"];
lQ = length(Qualities);
%%
for i = 1:ldays
    day_ = days(i);
    [yyyy,mm,dd] = datevec(day_);

    %%
    yyyyStr = sprintf("%04d",yyyy);
    mmStr = sprintf("%02d",mm);
    ddStr = sprintf("%02d",dd);

    %%
    for j = 1:lQ
        Q_ = Qualities(j);
        if strcmp(Q_,"NRTI")
            dirName = sprintf("S5P_DLR_%s_01_L3_SO2_%04d%02d%02d",Q_,yyyy,mm,dd);
            fName = sprintf("%s.nc",dirName);
            status = download_tropomi(dirName,fName,yyyyStr,mmStr,ddStr);
            if ~status
                break;
            end
            dirName = sprintf("S5P_DLR_%s_01_040201_L3_SO2_%04d%02d%02d",Q_,yyyy,mm,dd);
            fName = sprintf("%s.nc",dirName);
            status = download_tropomi(dirName,fName,yyyyStr,mmStr,ddStr);
            if ~status
                break;
            end
            dirName = sprintf("S5P_DLR_%s_01_040100_L3_SO2_%04d%02d%02d",Q_,yyyy,mm,dd);
            fName = sprintf("%s.nc",dirName);
            status = download_tropomi(dirName,fName,yyyyStr,mmStr,ddStr);
            if ~status
                break;
            end
        elseif strcmp(Q_,"RPRO")
            dirName = sprintf("S5P_DLR_%s_01_L3_SO2_%04d%02d%02d",Q_,yyyy,mm,dd);
            fName = sprintf("%s.nc",dirName);
            status = download_tropomi(dirName,fName,yyyyStr,mmStr,ddStr);
            if ~status
                break;
            end
        elseif strcmp(Q_,"OFFL")
            dirName = sprintf("S5P_DLR_%s_01_L3_SO2_%04d%02d%02d",Q_,yyyy,mm,dd);
            fName = sprintf("%s.nc",dirName);
            status = download_tropomi(dirName,fName,yyyyStr,mmStr,ddStr);
            if ~status
                break;
            end
        end
    end
end
end

function status = download_tropomi(dirName,fileNames,yyyyStr,mmStr,ddStr)
status = 1;
lFiles = length(fileNames);

for k = 1:lFiles
    file_ = fileNames(k);
    cmd = strcat("curl -f --remote-name 'https://download.geoservice.dlr.de/S5P_TROPOMI/files/L3/",...
        yyyyStr,"/",mmStr,"/",ddStr,"/",dirName,'/',file_,"'");
    fprintf("%s\n",cmd);

    %%
    status = unix(cmd);
    fprintf("\n");
    fprintf("status: %d\n",status);
    fprintf("\n");
    if ~status
        return;
    end
end
end