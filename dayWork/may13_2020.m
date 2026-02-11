clear; close all; clc;

baseDir = "~/masa/old/research/now/focmechs/";
cd(baseDir);
load SierraNegraPolaritiesThirdSet.mat;
E = table2struct(E);
E = E(13);
newDepths = E.depth; %(-1.1:0.1:3.9)';

%%
for i = 1:length(newDepths)
    E_ = E;
    E_.depth = newDepths(i);

    [strike_,dip_,rake_] = SC2HASH(E_,false,false);
    strike(i,1) = strike_;
    dip(i,1) = dip_;
    rake(i,1) = rake_;

    [strike2_,dip2_,rake2_] = auxplane(strike_,dip_,rake_);
    strike2(i,1) = strike2_;
    dip2(i,1) = dip2_;
    rake2(i,1) = rake2_;

    %%
    mt = sdr2mt(strike_,dip_,rake_);
    ax = plotEvents(E_,true);
    hold(ax,'on');
    c = colorbar;

    ax(2) = axes;
    focalmech(mt,E_.lon,E_.lat,0.01);
    %focalmech(ax(2),mt,E_.lon,E_.lat,0.01);

    ax(2).Visible = 'off';
    axis(ax(2),'equal');
    linkaxes(ax,'xy');
    axis(ax,[E_.lon - 0.1,E_.lon + 0.1,E_.lat - 0.1,E_.lat + 0.1]);

    id = E_.id;
    t = E_.t;
    mag = E_.mag;
    lat = E_.lat;
    lon = E_.lon;
    depth = E_.depth;

    cd(baseDir);
    cd(id)

    depth_ = 1e3*(newDepths(i)+1.1);
    depthStr = num2str(depth_);

    if depth_ < 100
        depthStr = ['000' depthStr];
    elseif depth_ < 1000
        depthStr = ['0' depthStr];
    end

    fname = strcat(id,"_",num2str(depthStr),"_TeppHashFocalMechanism");
    title(ax(1),{...
        strcat(datestr(t),", mag: ",num2str(mag),...
        ", lat: ",num2str(lat),", lon: ",num2str(lon),...
        ", depth: ",num2str(depth)); ...

        strcat("NP1: ",num2str(strike(i)),...
        ", ",num2str(dip(i)),", ",num2str(rake(i)),...
        "; NP2: ",num2str(round(strike2(i))),...
        ", ",num2str(round(dip2(i))),", ",num2str(round(rake2(i))))},...
        'interpreter','latex');
    disp(fname);
    print('-djpeg',fname);
end