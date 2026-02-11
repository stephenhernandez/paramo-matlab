function [t,lat,lon,depth,mag,...
    timerr,laterr,lonerr,deptherr,magerr,...
    rms,azgap,usedPhases,nmag,...
    id,methodID,earthModel,evMode,evStatus,agencyID,magtype,magmethod,evDescription,...
    creationTime,authorID,evType] = readSourceParameters(event_file)

%%
fid = fopen(event_file);
n = 0; %number of lines read
while n < 1
    buff = fgetl(fid);
    if buff ~= -1
        n = n+1;
        tmp = textscan(buff,'%s');
        tmp = tmp{1};
        ltmp = length(tmp);
        if ltmp ~= 36
            disp('not all fields were read. something is wrong. double check. sorry.');
        else
            data = NaN(19,1);
            for i = 1:16
                data(i) = str2double(tmp{i});
            end
            t = datetime(data(1:6)');
            timerr = data(7);
            lat = data(8);
            laterr = data(9);
            lon = data(10);
            lonerr = data(11);
            depth = data(12);
            deptherr = data(13);
            rms = data(14);
            azgap = data(15);
            usedPhases = data(16);
            id = tmp{17};
            methodID = tmp{18};
            earthModel = tmp{19};
            evMode = tmp{20};
            evStatus = tmp{21};
            agencyID = tmp{22};
            for i = 23:25
                data(i-6) = str2double(tmp{i});
            end
            mag = data(17);
            magerr = data(18);
            nmag = data(19);
            magtype = tmp{26};
            magmethod = tmp{27};
            evDescription = tmp{28};
            creationTime = datetime(str2double(tmp{29}),str2double(tmp{30}),str2double(tmp{31}),...
                str2double(tmp{32}),str2double(tmp{33}),str2double(tmp{34}));
            authorID = tmp{35};
            evType = tmp{36};
        end
    end
end
fclose(fid);
