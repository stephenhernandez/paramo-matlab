function [t,lambda,angs,vangs,rectilinearity,planarity] = polarizeLoop(S,twin,lfc,hfc,diffFlag)
if nargin < 3; lfc = -inf; end
if nargin < 4; hfc = -inf; end
if nargin < 5; diffFlag = false; end

if diffFlag
    S = differentiateWaveforms(S);
end
S = syncWaveforms(S);
Fs = round(1/S(1).delta);
data = double(pull(S));

npts=size(data,1);
data = data(1:npts,:);
data = detrend(data);
e = data(:,1);
n = data(:,2);
z = data(:,3);
t = getTimeVec(S(1));

npoles = 4;
zeroPhaseFlag = false;

%% filter the data
if any(isfinite([lfc hfc]))
    disp('filtering data')
    disp([lfc hfc])
    df = zpkFilter(data,lfc,hfc,Fs,npoles,zeroPhaseFlag);
    ef = df(:,1);
    nf = df(:,2);
    zf = df(:,3);
    clear df;
else
    ef = e;
    nf = n;
    zf = z;
end

%% create 6 streams
e2 = ef.^2;
n2 = nf.^2;
z2 = zf.^2;
en = ef.*nf;
ez = ef.*zf;
nz = nf.*zf;
data = [e2 n2 z2 en ez nz];

%% get variance-covariance coeffs (assuming zero-mean)
disp('getting coefficients');
tic;
winlen = twin*Fs;
boxcar = ones(winlen,1);
data = fftfilt(boxcar,data);
toc;

%% get the job done
t2 = t;
Ndata = length(t2);
data = data(1:Ndata,:);
angs = NaN(Ndata,3);
vangs = angs;
lambda = angs;
c = NaN(3);

tic;
disp('performing time loop')
for i = 1:Ndata
    %populate variance-covariance matrix
    c(1,1) = data(i,1); c(1,2) = data(i,4); c(1,3) = data(i,5);
    c(2,1) = data(i,4); c(2,2) = data(i,2); c(2,3) = data(i,6);
    c(3,1) = data(i,5); c(3,2) = data(i,6); c(3,3) = data(i,3);

    [v1,d1]=eig(c);      	        % eigenvalue/eigenvectors
    [eigVectors,d]=order3(v1,d1);  	% order the eigenvalues and eigenvectors

    if i == 4370
    disp(eigVectors(:,1))
    rssq(eigVectors(:,1))
    rssq(eigVectors(1:2,1))
    end
    %disp(rssq(eigVectors))
    %azimuths
    angs(i,1) = 90 - atan2d(eigVectors(2,1),eigVectors(1,1)); % azimuth for first  eigenvalue
    angs(i,2) = 90 - atan2d(eigVectors(2,2),eigVectors(1,2)); % azimuth for second eigenvalue
    angs(i,3) = 90 - atan2d(eigVectors(3,3),eigVectors(1,3)); % azimuth for third  eigenvalue

    %eigenvalues
    lambda(i,1) = d(1);
    lambda(i,2) = d(2);
    lambda(i,3) = d(3);

    %angle from the vertical
    vangs(i,1) = 90 - asind(eigVectors(3,1)); % 180/pi;
    vangs(i,2) = 90 - asind(abs(eigVectors(3,2))); %* 180/pi;
    vangs(i,3) = 90 - asind(abs(eigVectors(3,3))); %* 180/pi;
end
disp('done with time loop');
toc;

%%
aI = angs >= 360;
angs(aI) = angs(aI) - 360;
aI = angs < 0;
angs(aI) = angs(aI) + 360;
angs = angs - 180;

jurkFlag = ~true;
if jurkFlag
    rectilinearity = 1-((lambda(:,2) + lambda(:,3))./(2*lambda(:,1)));  % jurkevics 88
    planarity = 1-((2*lambda(:,3))./(lambda(:,1)+lambda(:,2)));         % jurkevics 88
else
    rectilinearity = 1-((lambda(:,2) + lambda(:,3))./lambda(:,1));      % vidale 86
    planarity = 1-(lambda(:,3)./lambda(:,2));                           % vidale 86
end

%% figure 1
plotFlag= true;
if plotFlag
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(4,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    ha3(1) = nexttile;
    lambda = abs(lambda);
    plot(t2,lambda);
    set(gca,'YScale','log')
    axis tight
    grid on;
    tstr = 'raw eigenvalues';
    title(tstr);

    ha3(2) = nexttile;
    plot(t2,lambda./lambda(:,1)); %hold on; plot(t2,de(:,2)./de(:,1)); plot(t2,de(:,3)./de(:,1));
    axis tight
    grid on;
    title('relative eigenvalues');

    ha3(3) = nexttile;
    plot(t,ef);
    grid on;
    title('EW Comp');

    ha3(4) = nexttile;
    plot(t2,angs(:,1),'.','markersize',18);
    axis tight
    grid on;
    title('azimuth ');

    ha3(5) = nexttile;
    plot(t,nf);
    grid on;
    title('NS Comp');

    ha3(6) = nexttile;
    plot(t2,vangs(:,1));
    grid on;
    title('angle from vertical'); %incidence angle ');
    %ylim([0 90]);

    ha3(7) = nexttile;
    plot(t,zf);
    grid on;
    title('Z comp');

    ha3(8) = nexttile;
    plot(t2,rectilinearity);
    grid on;
    hold on;
    plot(t2,planarity); zoom on;
    grid on;
    if jurkFlag
        title("Jurkevics ('88)");
    else
        title("Vidale ('86)");
    end
    legend('lin.','plan.','location','best');
    linkaxes(ha3,'x');
    sgtitle(S(1).kstnm);

    % %% figure 2, plot the raw data
    % figure('units','normalized','outerposition',[0 0 1 1]);
    % ha(1) = subplot(311);
    % plot(t,e);
    % grid on;
    % title('EW Comp');
    % 
    % ha(2) = subplot(312);
    % plot(t,n);
    % grid on;
    % title('NS Comp');
    % 
    % ha(3) = subplot(313);
    % plot(t,z);
    % grid on;
    % title('Z comp');
    % linkaxes(ha,'x');
    % sgtitle(S(1).kstnm);
    % zoom on;
end
