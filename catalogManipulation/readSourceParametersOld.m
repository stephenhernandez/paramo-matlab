function [origyyyy,origmm,origdd,orighh,origmmm,origsec,origlat,origlaterr,...
    origlon,origlonerr,origdepth,origrms,origgap,id] = readSourceParametersOld(event_file)

fid = fopen(event_file);
eofstat = feof(fid);
idFlag = false;

%%
while ~eofstat
    buff = fgetl(fid);
    if buff ~= -1
        tmp = textscan(buff,'%s');
        tmp = tmp{1}; % read first line
        if strcmp(tmp{1},'Public')
            if ~idFlag
                C = textscan(buff, '%s %s %s');
                id = C{3};
                id = id{1};
                idFlag = true;
            end
        elseif strcmp(tmp{1},'Date')
            origyyyy = sscanf(buff(28:31),'%f');
            origmm = sscanf(buff(33:34),'%f');
            origdd = sscanf(buff(36:37),'%f');
        elseif strcmp(tmp{1},'Time')
            orighh = sscanf(buff(28:29),'%f');
            origmmm = sscanf(buff(31:32),'%f');
            origsec = sscanf(buff(34:39),'%f');
        elseif strcmp(tmp{1},'Latitude')
            C = textscan(buff,'%s %f %s %s %f %s');
            origlat = C{2};
            err_ = C{5};
            if isempty(err_)
                origlaterr = NaN;
            else
                origlaterr = err_;
            end
        elseif strcmp(tmp{1},'Longitude')
            C = textscan(buff,'%s %f %s %s %f %s');
            origlon = C{2};
            err_ = C{5};
            if isempty(err_)
                origlonerr = NaN;
            else
                origlonerr = err_;
            end
        elseif strcmp(tmp{1},'Depth')
            C = textscan(buff,'%s %f %s %s');
            origdepth = C{2};
        elseif strcmp(tmp{1},'Residual')
            C = textscan(buff,'%s %s %f %s');
            origrms = C{3};
        elseif strcmp(tmp{1},'Azimuthal')
            C = textscan(buff,'%s %s %f %s');
            origgap = C{3};
        end
    end
    eofstat = feof(fid);
end
fclose(fid);
