%function [observable,t0,duration] = readObservable2(stnm,enz,lfc,hfc,npoles,filtFlag,flag3,fs,t0,baz,secDur)
function [observable,t0,duration] = readObservable2(stnm,enz,lfc,hfc,npoles,filtFlag,flag3,fs,t0,az,maxDur)
delta = 1/fs;

if filtFlag
    Hd = zpkOperator(lfc,hfc,fs,npoles);
    enz = detrend(filter(Hd,enz));
end

if flag3
    oz=enz(:,3);
    on=enz(:,2);
    oe=enz(:,1);

    % Cut the data
    t = ((0:length(oz)-1)*delta)';
    endTry = t0 + maxDur;
    tI = t >= t0 & t <= endTry;

    oz = detrend(oz(tI));
    on = detrend(on(tI));
    oe = detrend(oe(tI));
    if ~az
        [tobs, robs] = rotate2d(oe,on,az);
    end
else
    oz=enz; %(:,3);

    % Cut the data
    t = ((0:length(oz)-1)*delta)';
    endTry = t0 + maxDur;
    tI = t >= t0 & t <= endTry;

    oz = detrend(oz(tI));
end

%%
if flag3
    if isempty(oz) || isempty(robs) || isempty(tobs)
        fprintf('Skipping: %s\n',stnm);
        observable = [];
        duration=length(observable)/3;
    else
        %observable = [oz;robs;tobs];
        %observable = [tobs;robs;oz];
        if strcmp(stnm,'ACHA1')
            observable = [zeros(size(tobs));robs;oz];
        elseif strcmp(stnm,'POPE1')
            observable = [zeros(size(tobs));robs;oz];
        else
            observable = detrend([tobs;robs;oz]);
        end
        duration=length(observable)/3;
        %duration=length(observable)/2;
    end
else
    observable = oz;
    duration=length(observable);
    %     if isempty(oz) || isempty(robs) || isempty(tobs)
    %         disp(['Skip ',stnm])
    %         observable = [];
    %         duration=2*length(observable)/3;
    %     else
    %         %observable = [oz;robs;tobs];
    %         %observable = [tobs;robs;oz];
    %         if strcmp(stnm,'POPE1')
    %             observable = [zeros(size(tobs));robs];
    %         elseif strcmp(stnm,'ACHA1')
    %             observable = [zeros(size(tobs));robs];
    %         else
    %             %observable = [tobs;robs]; duration=length(observable)/3;
    %             %observable = [robs;oz]; duration=length(observable)/2;
    %             observable = oz;
    %             duration=length(observable);
    %         end
    %     end
    %     %
    %     %         %just the vertical
    %     %         if isempty(oz)
    %     %             observable = [];
    %     %             duration=length(observable);
    %     %         else
    %     %             %observable = oz;
    %     %             observable = robs;
    %     %             %observable = tobs;
    %     %             duration=length(observable);
    %     %         end
end
%disp([stnm ' ' num2str(baz) ' ' num2str(duration)])
