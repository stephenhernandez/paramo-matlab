function updateGlobalCatalog(lastDay)
%cd ~/igdata/
load('~/igdata/globalCatalog','t','eqlat','eqlon','eqdepth','eqmag','magType','status','code');

%%
td = datenum(t);
tomorrowDay = ceil(now);
if nargin < 1
    lastDay = floor(max(td))-2;
end


%%
lbc = LibComCat();
events = lbc.getEventData('starttime',lastDay,'endtime',tomorrowDay, ...
    'xmin',-180,'xmax',180,'ymin',-90,'ymax',90,'minmag',5,'maxmag',10);
U = libcomcat2struct(events);
t_ = pull(U,'t');
t = [t; t_];

eqlat_ = pull(U,'lat');
eqlat = [eqlat; eqlat_];

eqlon_ = pull(U,'lon');
eqlon = [eqlon; eqlon_];

eqmag_ = pull(U,'mag');
eqmag = [eqmag; eqmag_];

eqdepth_ = pull(U,'depth');
eqdepth = [eqdepth; eqdepth_];

magType_ = string(pull(U,'magType'));
magType = [magType; magType_];

status_ = string(pull(U,'status'));
status = [status; status_];

code_ = string(pull(U,'code')); 
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
clearvars -except t eqlat eqlon eqdepth eqmag magType status code
save('~/igdata/globalCatalog','t','eqlat','eqlon','eqdepth','eqmag','magType','status','code');
