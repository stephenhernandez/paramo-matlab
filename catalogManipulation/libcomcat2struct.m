function U = libcomcat2struct(events)
levents = length(events);
U = populateUSGSCatStructure(levents);

for i = 1:levents
    properties = events(i).properties;
    coords = events(i).geometry.coordinates;
    U(i).mag = properties.mag;
    U(i).lat = coords(2);
    U(i).lon = coords(1);
    U(i).depth = coords(3);
    U(i).place = properties.place;
    U(i).t = datetime(1970,01,01) + milliseconds(properties.time); %datetime(properties.time/1000,'ConvertFrom','posixtime');
    U(i).updated = datetime(1970,01,01) + milliseconds(properties.updated); %datetime(properties.updated/1000,'ConvertFrom','posixtime');
    U(i).tz = properties.tz;
    U(i).url = properties.url;
    U(i).detail = properties.detail;
    U(i).felt = properties.felt;
    U(i).cdi = properties.cdi;
    U(i).mmi = properties.mmi;
    U(i).alert = properties.alert;
    U(i).status = properties.status;
    U(i).tsunami = properties.tsunami;
    U(i).sig = properties.sig;
    U(i).net = properties.net;
    U(i).code = properties.code;
    U(i).ids = properties.ids;
    U(i).sources = properties.sources;
    U(i).types = properties.types;
    U(i).nst = properties.nst;
    U(i).dmin = properties.dmin;
    U(i).rms = properties.rms;
    U(i).gap = properties.gap;
    magType = properties.magType;
    if isempty(magType)
        magType = '';
    end
    U(i).magType = magType;
    U(i).type = properties.type;
    U(i).title = properties.title;
end