function plotRecordSection(S,varargin)

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Aug 13, 2019

%%
nVarargin = length(varargin);
functionDefaults = {...
    -inf,...                        % lfc
    -inf,...                        % hfc
    true,...                        % vFlag
    true,...                        % saveFlag
    [],...                          % default id
    'Z',...                         % cmp
    false};                         % diffFlag

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[lfc,hfc,vFlag,saveFlag,id,cmp,diffFlag] = deal(optsToUse{:});

%%
PWD = pwd;
if isempty(lfc)
    lfc = -inf;
end

if isempty(hfc)
    hfc = -inf;
end

%%
S = S(:,1);

%%
if diffFlag
    S = differentiateWaveforms(S);
end

%%
cornersfin = isfinite([lfc hfc]);
if any(cornersfin)
    npoles = 4;
    zeroPhaseFlag = false;
    tw = 0.02;
    S = filterWaveforms(S,lfc,hfc,npoles,tw,zeroPhaseFlag);
end

%%
figure('units','normalized','outerposition',[0 0 1 1]);
ax = newplot();
h = gcf;

if ~vFlag
    h.Visible = 'off';
end

%%
charchan = char(pull(S,'kcmpnm'));
cmpI = charchan(:,3) == cmp;
dist = pull(S,'dist');
cmpI = cmpI & isfinite(dist);
sumCmps = sum(cmpI);

%%
if isempty(id)
    id = char(S(1).evid);
    if isempty(id)
        id = datestr(S(1).ref,30);
    end
end

%%
eqlat = S(1).evla;
eqlon = S(1).evlo;
eqdepth = S(1).evdp;
eqmag = S(1).eqmag;

%%
fontSize = 18;
if sumCmps
    hold(ax,'on');
    %S = normalizeWaveforms(S,true,true);
    distz = dist(cmpI);
    ampFact = ceil(max(distz))*5e-6;
    S = S(cmpI);
    for j = 1:sumCmps
        ax.ColorOrderIndex = j;
        plot(ax,getTimeVec(S(j)),ampFact*S(j).d+distz(j),'LineWidth',1);
        text(ax,S(j).ref+S(j).e,distz(j),S(j).kstnm,'FontSize',fontSize);
    end
    ylabel(ax,'Distance [km.]');
    titleStr = {['lat: ',num2str(eqlat),', lon: ',num2str(eqlon),', depth: ',num2str(eqdepth)];...
        ['id: ',char(id),', mag: ',num2str(eqmag)]};
    title(ax,titleStr);
    
    if saveFlag
        cd ~/products/events/html/
        if ~exist(id,'dir')
            mkdir(id);
        end
        cd(id)
        print('-djpeg',strcat(id,'_recordSection'));
        close(h);
    else
        zoom(ax,'on');
    end
else
    fprintf(2,'no data found\n');
end
cd(PWD);
