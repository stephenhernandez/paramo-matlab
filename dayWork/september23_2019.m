clear; close all; clc;
%cd /home/shernandez/products/duplicates
%load ispt2016locatedtemplates.mat
cd ~/research/now/reve_silvi/
%load vent_e_infrasound_template
%load vent_e_9_infrasound_template
%load negative_polarity_template
% template = flipud(template);
% template = template(1:1500);
% template = flipud(template);

load vent_e_9_seismic_template
template = flipud(template);
template = decimate(template,2);
template = template(1:3000);
template = template/norm(template);

%%
template = flipud(template);
warning off signal:findpeaks:largeMinPeakHeight

%%
pre = 0;
post = 15;
days = datetime(2018,01,01):datetime(2018,12,31);
lDays = length(days);
pks = cell(lDays,1);
locs = pks;
mads = pks;
lfc = 1/6;
hfc = 1/3;
winlen = length(template);
matches = NaN(winlen,5000);
n = 0;
for i = 1:lDays
    disp(i)
    tStart = days(i);
    [yyyy_,mm_,ddd_] = datevec(tStart);
    S = loadWaveforms(tStart,1,"REVN","HHZ","EC");

    if any(isnat(pull(S,'ref')))
        continue
    else
        S = intWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc));
        S = syncWaveforms(S);
        S = resampleWaveforms(S,50);
        Fs = 1./S.delta;
        data = double(pull(S));
        
        winlen = size(template,1);
        lData = size(data,1);
        
        %%
        if lData >= winlen
            
            data2 = double(data).^2;
            box = ones(winlen,length(S));
            norms = fftfilt(box,data2);
            norms = sqrt(abs(norms));
            t = getTimeVec(S(1));
            
            ccnorm = fftfilt(double(template),data);
            snorms = norms > 1;
            ccnorm = snorms.*ccnorm./norms;
            
            stack = nanmean(ccnorm,2); mad_ = nanstd(stack);
            newThresh_ = min([0.7 15*mad_]);
            [pks_,locs_] = findpeaks(stack,'MinPeakHeight',newThresh_,'MinPeakDistance',winlen);
            locs_ = locs_-winlen+1;
            locsI = locs_ > 0;
            pks_ = pks_(locsI);
            locs_ = locs_(locsI);
            lpks = length(pks_);
            
            %%
            if lpks
                fprintf('Found %d event(s)\n',lpks);
                pks{i} = pks_;
                mads{i} = pks_/mad_;
                locs{i} = datenum(t(locs_));
                for ii = 1:lpks
                    n = n+1;
                    newData = data(locs_(ii):locs_(ii)+winlen-1,:);
                    lnewData = size(newData,1);
                    matches(1:lnewData,n) = newData;
                end
            end
            
        end
    end
end
matches = matches(:,1:n);

%%
pksOrig = cat(1,pks{:});
tnew = dn2dt(cat(1,locs{:}));
maxAmpOrig = rms(matches)';

badI = (pksOrig < 0.85);
tnew(badI) = [];
matchesnew = matches(:,~badI);
maxAmp = rms(matchesnew)';
pksOrig(badI) = [];

close all;
figure(1); hold on;
S = scatter(tnew,(1:length(tnew))',3*exp(log10(maxAmp)),pksOrig,'d','filled');
zoom on; colorbar; caxis([0.85 1])
S.MarkerFaceAlpha = 0.4;
S.MarkerEdgeColor = 'k';
xlim([datetime(2018,01,01) datetime(2018,12,31)]);

tdum = seconds((0:winlen-1)'/Fs);
figure(2); hold on;
plot(tdum,normalizeWaveforms(matchesnew));
hold on;
plot(tdum,normalizeWaveforms(nanmean(matchesnew,2)),'k','linewidth',4);
axis tight;
zoom on;

figure(3); 
hold on;
S = scatter(tnew,log10(maxAmp),8*exp(log10(maxAmp)),pksOrig,'d','filled');
c = colorbar; caxis([0.85 1]); ylim([1.5 3]);
S.MarkerFaceAlpha = 0.4;
S.MarkerEdgeColor = 'k';
zoom on;
ylabel('$log_{10}(A)$');
c.Label.String = 'Correlation Coefficient';
c.Label.Interpreter = 'latex';
xlim([datetime(2018,01,01) datetime(2018,12,31)]);

%%
% tI = [];
% tMaster = [];
% pksMaster = [];
% for i = 1:lT
%     tnew = dn2dt(cat(1,locs{:,i}));
%     tMaster = [tMaster; tnew];
%     pksMaster = [pksMaster; cat(1,pks{:,i})];
%     tI = [tI; i*ones(length(tnew),1)];
% end
% [tMaster,sI] = sort(tMaster);
% pksMaster = pksMaster(sI);
% tI = tI(sI);
%
%
% [t_,cc_,removeIndices] = removeRepeatedMatches(tMaster,pksMaster,30);
% for jj = 1:length(removeIndices)
%     tI(removeIndices{jj}) = [];
% end
%
% ccI = cc_ >= 0.65;
% t_ = t_(ccI);
% cc_ = cc_(ccI);
% tI = tI(ccI);
%
% close all
% figure(); S = scatter(t_,1:length(t_),[],cc_,'o','filled'); zoom on; colorbar;