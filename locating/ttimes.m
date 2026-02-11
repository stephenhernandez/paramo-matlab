% function [X,T] = ttimes(h0,wavetype,model)
% Calculate travel time to surface for a given take-off angle and source
% depth
%
% Heavily modified (by Stephen Hernandez) version of code from UW group

%%
function [X, T] = ttimes(h0,wavetype,model)
%Define Costa Rica default wavetype, and velocity mode

if nargin < 2
    wavetype = 's';
end

if nargin < 3
    model = [
        3.01 0.00 5.35
        3.44 5.00 6.12
        3.53 13.00 6.28
        3.63 16.00 6.46
        3.78 20.00 6.72
        3.94 25.00 7.01
        4.15 30.00 7.39
        4.24 35.00 7.55
        4.57 40.00 8.14
        4.64 65.00 8.26];
end


if strcmp(wavetype,'p')
    vm = model(:,3);
else
    vm = model(:,1);
end

vmax = max(vm);
vmin = min(vm);

h = model(:,2);
hmax = max(h);

takeoffAngles = flipud(90 - logspace(log10(0.01),log10(89.99),600)'); %[(0.05:0.05:89.85)'; 89.9; 89.95]; %downgoing rays
turncrit = true(size(takeoffAngles));

%Isolate layers above and below source depth
if h0 >= hmax
    takeoffAngles = [-flipud(takeoffAngles); 0]; % only up-going rays
else
    takeoffAngles = [takeoffAngles; -flipud(takeoffAngles)]; % down-going and up-going rays
    turncrit = [~turncrit; turncrit];
end
takeoffAngles = deg2rad(takeoffAngles);

hI = h0 > h;
htop = [h(hI); h0];
dhtop = diff(htop);
hbot = [h0; h(~hI)];
dhbot = diff(hbot);

vtop = [vm(hI); vm(sum(hI))];
dvtop = diff(vtop);
vbot = [vm(sum(hI)); vm(~hI)];
dvbot = diff(vbot);

% preallocation
X = zeros(length(takeoffAngles),1);
T = X;
lTOA = length(takeoffAngles);

%get job done
for i = 1:lTOA
    n = 1;
    v = vbot(1);
    theta = takeoffAngles(i);   %initial take-off angle
    p = sin(theta)/v;           %ray parameter, snells law

    if p > 0                    %Down-going rays
        dv = dvbot;
        dh = dhbot;
        while v < vmax && n <= length(dh)
            iold = theta;
            vold = v;
            X(i) = X(i)+abs(dh(n)*tan(theta));

            T(i) = T(i)+abs(dh(n)/(v*cos(theta)));
            v = v+dv(n);
            theta = asin(v*p);

            if ~isreal(theta)
                theta = -iold ;
                v = vold;
                p = -p;
                dh = [dh(1:n); flipud(dh(1:n)); flipud(dhtop)];
                dv = [dv(1:n); -flipud(dv(1:n-1)); -flipud(dvtop); 0];
                turncrit(i) = true;
                %disp(turncrit(ii))
            end
            n = n+1;
        end
    else                        %Up-going rays
        dv = [-flipud(dvtop(1:end-1)); 0];
        dh = flipud(dhtop);
        while v >= vmin && n <= length(dh)
            X(i) = X(i)+abs(dh(n)*tan(theta));
            T(i) = T(i)+abs(dh(n)/(v*cos(theta)));
            v = v+dv(n);
            theta = asin(v*p);
            n = n+1;
        end
    end
end

X = X(turncrit);
T = T(turncrit);
takeoffAngles = takeoffAngles(turncrit);
%figure; plot(X,T,'ko-');
%find up-going rays
Xold = X;
Told = T;
%Aold = angles;
ugr = takeoffAngles < 0;
data = [Xold(ugr),Told(ugr)];%,Aold(ugr)];
data = sortrows(data);
Xup = data(:,1); Tup = data(:,2); %Aup = data(:,3);
XI = Xup<500;
X = [Xup(XI); 500];
Tup = interp1(Xup,Tup,X);

%find down-going rays
dgr = takeoffAngles >= 0;
data = [Xold(dgr),Told(dgr)];
data = sortrows(data);
Xdown = data(:,1);
Tdown = interp1(Xdown,data(:,2),X);
T = min([Tdown';Tup'])';
