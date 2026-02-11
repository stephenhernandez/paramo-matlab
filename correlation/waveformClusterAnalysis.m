function [linear_stack,singletons] = waveformClusterAnalysis(indiv_events,tabs,thresh,Fs,plotFlag)
if nargin < 3; thresh = 0.7; end
if nargin < 4; Fs = 100; end
if nargin < 5; plotFlag = true; end

%load('wolf_manually_picked');
maxAmpRMS = rms(indiv_events);
[maxccp,plags,maxccn,nlags] = docc(indiv_events,true);
nLarger = maxccn > maxccp;
normcc_down = maxccp;

tic;
normcc_down = squareform(normcc_down);
toc;
disp('done getting correlation matrix');
nSearch = length(normcc_down);

if plotFlag
    figh(6) = figure('units','normalized','outerposition',[0 0 1 1]);
    h = imagesc(normcc_down);
    title('Similarity Matrix');
    ylabel('Event Number');
    xlabel('Event Number');
    set(h, 'alphadata', normcc_down >= thresh);
    axis square;
    colorbar;
    caxis([thresh 1]);
end

disp('counting number of repeats')
nMatches = zeros(nSearch,1);
n_flipped = 0;
flipped_indices = cell(nSearch,1);
nLSquare = squareform(nLarger);
ccDiff = squareform(maxccn-maxccp);
diffThresh = 0.15;
for i = 1:nSearch
    normcc_ = normcc_down(:,i);
    nLSquare_ = nLSquare(:,i);
    nI = normcc_ >= thresh;
    LOCS_ = find(nI);
    nMatches(i) = length(LOCS_);
    if ~isempty(LOCS_) % if found events other than yourself, then ...?
        disp(['Event ',num2str(i),'; repeats: ',num2str(nMatches(i)),', ',datestr(tabs(i))])
        LOCS2_ = ccDiff(:,i)>= diffThresh;
        posThresh1Thresh2 = nLSquare_&nI&LOCS2_;
        sumnLSquare_ = sum(posThresh1Thresh2);
        if sumnLSquare_
            n_flipped = n_flipped + sumnLSquare_;
            flipped_indices{i} = find(posThresh1Thresh2);
            disp(['    Number of possible flips: ',num2str(sumnLSquare_)])
            disp(['    Total number of possible flips so far: ',num2str(n_flipped)])
            disp('    Index(es) of possible flips: ');
            disp(flipped_indices{i})
        end
    end
end
flipped_indices{nSearch+1} = unique(sort(cat(1,flipped_indices{1:nSearch}))); %store the data
disp('Display index(es) of all possible flips...');
disp(flipped_indices{nSearch+1})
clear normcc_ pol_ nI LOCS_ master posThresh1Thresh2 LOCS2_ ccDiff

if plotFlag
    markersize = 15;
    figh(7) = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(211)
    plot(nMatches,'c.','markersize',markersize);
    xlabel('Event Number');
    ylabel('Number of Repetitions');
    subplot(212)
    plot(tabs,nMatches,'r.','markersize',markersize);
    ylabel('Number of Repetitions');
    ax = gca;
    ax.XTickLabelRotation = 45;
end

normcc_down = 1-normcc_down; %turn similarity matrix into distance matrix
dI = normcc_down > 0.99999;
normcc_down(dI) = 0;
normcc_down = squareform(normcc_down); %convert distance matrix to distance vector
clear dI
method = 'average';
[family,l_uniq_indices] = returnClusts(normcc_down,thresh,method);
%ClustTot = length(l_uniq_indices);
%mostMembers = l_uniq_indices(1);
singletonI = l_uniq_indices < 2;
nSingletons = sum(singletonI);
nS = 0;
for i = 1:length(l_uniq_indices)
    if l_uniq_indices(i) == 1
        nS = nS+1;
        family1 = family{i};
        singletons(:,nS) = indiv_events(:,family1);
    end
end

l_families = 0;
maxNFam = 750;
nFam = sum(~singletonI);
linear_stack = indiv_events(:,1:nFam);

maxccp = squareform(maxccp);
maxccn = squareform(maxccn);
plags = squareform(plags);
nlags = squareform(nlags);

%%
for i = 1:nFam
    family1 = family{i};
    lfam = length(family1);
    indiv_tmp = indiv_events(:,family1);
    disp(['Family ',num2str(i),': ',num2str(lfam)]);
    maxccptmp = maxccp(:,family1);
    maxccptmp = maxccptmp(family1,:);
    maxccptmp = squareform(maxccptmp)'; %column vector
    maxccntmp = maxccn(:,family1);
    maxccntmp = maxccntmp(family1,:);
    maxccntmp = squareform(maxccntmp)'; %column vector
    plagstmp = plags(:,family1);
    plagstmp = plagstmp(family1,:);
    plagstmp = squareform(plagstmp)';   %column vector
    nlagstmp = nlags(:,family1);
    nlagstmp = nlagstmp(family1,:);
    nlagstmp = squareform(nlagstmp)';   %column vector
    ccData = [maxccptmp plagstmp maxccntmp nlagstmp];
    disp(['ccData: ',num2str(size(ccData))])
    shifted_data = apply_vdcc(indiv_tmp,ccData);
    linear_stack(:,i) = normalize_traces(mean(shifted_data,2));
    if plotFlag
        if lfam >= 3
            l_families = l_families + 1;
            fig = figure(1000+i);
            fig.Units = 'normalized';
            fig.OuterPosition = [0 0 1 1];
            linear_stack(:,i) = plot_family_and_amps(shifted_data,tabs(family1),maxAmpRMS(family1),1:lfam,5,Fs,-inf,-inf); %,false,[refTime refTime+length(days)],[plotMinAmp plotMaxAmp]);
        end
    end
end