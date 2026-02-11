cd ~/masa/old/research/now/focmechs/
clear
load SierraNegraPolaritiesThirdSet.mat;
%load SierraNegraPolaritiesSecondSet.mat;

E = table2struct(E);
E = E(13);
idOrig = E.id;
close all; clc;
newDepths = (0:0.2:5)';
for i = 1:length(newDepths)
    close all;

    depth_ = 1e3*newDepths(i);
    if depth_ < 10
        depthStr = ['000',num2str(depth_)];
    elseif depth_ < 100
        depthStr = ['00',num2str(depth_)];
    elseif depth_ < 1000
        depthStr = ['0',num2str(depth_)];
    else
        depthStr = num2str(depth_);
    end
    id = [idOrig,'_',depthStr,'m'];
    E.id = id;
    E.depth = newDepths(i);

    [strike_,dip_,rake_,conf0_] = SC2HASH(E,true,true);
    strike(i,1) = strike_;
    dip(i,1) = dip_;
    rake(i,1) = rake_;
    conf0(i,1) = conf0_;
    [strike2_,dip2_,rake2_] = auxplane(strike_,dip_,rake_);
    strike2(i,1) = strike2_;
    dip2(i,1) = dip2_;
    rake2(i,1) = rake2_;
end