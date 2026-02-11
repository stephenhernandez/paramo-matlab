%fernandina foc mechs
clear; close all; clc;
cd ~/research/now/fernandina/
hypofile = 'fer_fms.txt';
unix(sprintf('cat %s | grep 0NLL > tmp_file.txt ',hypofile));
T = readtable('tmp_file.txt');
lT = size(T,1);
t = char(string(T.Var1));
grepString = string(t(:,1:12));
yyyy = str2double(string(t(:,1:4)));
mm = str2double(string(t(:,5:6)));
dd = str2double(string(t(:,7:8)));
HH = str2double(string(t(:,9:10)));
MM = str2double(string(t(:,11:12)));
SS = zeros(lT);
centiSS = SS;
lat = NaN(lT);
lon = lat;

for i = 1:lT
    sec_ = str2double(t(i,13:14));
    centisec_ = str2double(t(i,15:16));
    if ~isfinite(sec_)
        sec_ = 0;
    end
    sec_ = sec_ + centisec_/100;
    SS(i) = sec_;
end

%%
t = datetime([yyyy mm dd HH MM zeros(length(HH),1)]);
for i = 1:lT
    grepString_ = grepString(i);
    unix(sprintf('cat %s | grep -v "0NLL" | grep %s > polarity1.txt',hypofile,grepString_));
    P1 = readtable('polarity1.txt');
end