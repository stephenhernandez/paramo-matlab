function updateGlobalCatalog2(startTime,endTime)
%cd ~/igdata/
load('~/igdata/globalCatalog2','t','eqlat','eqlon','eqdepth','eqmag','magType','status','code');

%%
lbc = LibComCat();
events = lbc.getEventData('starttime',startTime,'endtime',endTime, ...
    'xmin',-180,'xmax',180,'ymin',-90,'ymax',90,'minmag',4,'maxmag',10);
U = libcomcat2struct(events);

%%
try
    lU = length(U);
    t_ = pull(U,'t');
    eqlat_ = pull(U,'lat');
    eqlon_ = pull(U,'lon');
    eqmag_ = pull(U,'mag');
    eqdepth_ = pull(U,'depth');
    magType_ = string(pull(U,'magType'));
    status_ = string(pull(U,'status'));
    code_ = string(pull(U,'code'));
catch
    return;
end

%%
if lU == 1 && (isnat(t_) || isnan(eqmag_))
    fprintf('no data found for day: %s\n',datestr(startTime));
    return;
end

t = [t; t_];
eqlat = [eqlat; eqlat_];
eqlon = [eqlon; eqlon_];
eqmag = [eqmag; eqmag_];
eqdepth = [eqdepth; eqdepth_];
magType = [magType; magType_];
status = [status; status_];
code = [code; code_];

%%
[t,ia] = unique(t);
eqlat = eqlat(ia);
eqlon = eqlon(ia);
eqdepth = eqdepth(ia);
eqmag = eqmag(ia);
magType = magType(ia);
status = status(ia);
code = code(ia);

ia = isfinite(t);
t = t(ia);
eqlat = eqlat(ia);
eqlon = eqlon(ia);
eqdepth = eqdepth(ia);
eqmag = eqmag(ia);
magType = magType(ia);
status = status(ia);
code = code(ia);

%%
fprintf('data length: <strong>%d</strong>\n',length(t));

%%
clearvars -except t eqlat eqlon eqdepth eqmag magType status code
save('~/igdata/globalCatalog2','t','eqlat','eqlon','eqdepth','eqmag','magType','status','code');
