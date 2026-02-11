function plot_waveform(synth,observable,duration_,dist,stnms,af,zf,fs,flag3,keepI)
if nargin < 6
    af = 1e2;
end

if nargin < 7
    zf = 2;
end

if nargin < 8
    fs = 2;
end
LineWidth = 3;

%%
scal = 1/norm(observable);
observable = scal*observable;
synth = scal*synth;

miny = round(min(dist))-20;
maxy = round(max(dist))+20;
miny = max([miny 0]);
maxdur = max(duration_/fs)+10;

figure('units','normalized','outerposition',[0 0 1 1]);
ns = length(duration_);
if ~flag3
    disp([duration_(1),ns]);
    e_n_z = reshape(observable,[duration_(1),ns]);
    s_e_n_z = reshape(synth,[duration_(1),ns]);
    hold on;
    ax = gca;
    for i = 1:ns
        tdum = (0:duration_(i)-1)/fs;
        
        %get components
        z = e_n_z(:,i);
        z_s = s_e_n_z(:,i);
        
        %plot observables
        %prs = 
        plot(tdum,zf*af*z + dist(i),'linewidth',LineWidth);  zoom on; grid on;%ylim([0 12])
        ax.ColorOrderIndex = 2;

        % plot synthetics
        plot(tdum,zf*af*z_s + dist(i),'linewidth',LineWidth); zoom on; grid on;
        
        % tighten axes
        ylim([miny maxy]);
        xlim([0 maxdur]);
        
        % meaningful labels
        xlabel('Time [sec.]')
        ylabel('Distance [km.]')
        title({'Vertical Component',[num2str(zf),'X']})
        text(maxdur-5,dist(i),stnms{i},'FontSize',20);
        ax.ColorOrderIndex = 1;
    end
else
    %disp([duration_(1),ns*3]);
    e_n_z = reshape(observable,[duration_(1),ns*3]);
    s_e_n_z = reshape(synth,[duration_(1),ns*3]);
    for i = 1:ns
        tdum = (0:duration_(i)-1)/fs;
        
        %get components
        e = e_n_z(:,1+(i-1)*3);
        e_s = s_e_n_z(:,1+(i-1)*3);
        n = e_n_z(:,2+(i-1)*3);
        n_s = s_e_n_z(:,2+(i-1)*3);
        z = e_n_z(:,3+(i-1)*3);
        z_s = s_e_n_z(:,3+(i-1)*3);
        
        ha(1) = subplot(131); hold on;
        plot(tdum,af*e + dist(i),'linewidth',LineWidth);  zoom on; grid on;%ylim([0 12])
        title({'Tangential Component'});
        ha(1).ColorOrderIndex = 2;

        ha(2) = subplot(132); hold on
        plot(tdum,af*n + dist(i),'linewidth',LineWidth);  zoom on; grid on;%ylim([0 12])
        title({'Radial Component'});
        ha(2).ColorOrderIndex = 2;

        ha(3) = subplot(133); hold on;
        plot(tdum,zf*af*z + dist(i),'linewidth',LineWidth);  zoom on; grid on;%ylim([0 12])
        ha(3).ColorOrderIndex = 2;
        %ha(3).ColorOrderIndex

        ha(1) = subplot(131); hold on;
        plot(tdum,af*e_s + dist(i),'linewidth',LineWidth);  zoom on; grid on;ylim([miny maxy]);
        xlim([0 maxdur]);
        xlabel('Time [sec.]');
        ylabel('Distance [km.]');

        ha(2) = subplot(132); hold on;
        plot(tdum,af*n_s + dist(i),'linewidth',LineWidth);  zoom on; grid on;ylim([miny maxy]);
        xlim([0 maxdur]);
        xlabel('Time [sec.]');

        ha(3) = subplot(133); hold on;
        prs = plot(tdum,zf*af*z_s + dist(i),'linewidth',LineWidth);  zoom on; grid on; ylim([miny maxy]);
        %ha(3).ColorOrderIndex = 2;
        xlim([0 maxdur]);
        xlabel('Time [sec.]')
        title({'Vertical Component',[num2str(zf),'X']})
        text(maxdur-5,dist(i),stnms{i},'FontSize',14);

        ha(1).ColorOrderIndex = 1;
        ha(2).ColorOrderIndex = 1;
        ha(3).ColorOrderIndex = 1;
    end
    linkaxes(ha,'xy');
end

%% plot zooms
af = 1;
zf = 1;
if flag3
    fig = figure;%('units','normalized','outerposition',[0 0 1 1]);
    figNum = fig.Number;
    close(figNum);
    disp([duration_(1),ns*3]);
    e_n_z = reshape(observable,[duration_(1),ns*3]);
    s_e_n_z = reshape(synth,[duration_(1),ns*3]);
    
    figNum_ = unique(ceil((1:ns)/3));
    for i = 1:length(figNum_)
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    
    for i = 1:ns
        %disp(i)
        figNum_ = ceil(i/3);
        rowNum_ = 3+i-figNum_*3;
        figure(figNum + figNum_ -1);
        %figure('units','normalized','outerposition',[0 0 1 1]);
        %figNum + figNum_ -1
        hold on;
        tdum = (0:duration_(i)-1)/fs;
        
        %get components
        e = e_n_z(:,1+(i-1)*3);
        e_s = s_e_n_z(:,1+(i-1)*3);
        n = e_n_z(:,2+(i-1)*3);
        n_s = s_e_n_z(:,2+(i-1)*3);
        z = e_n_z(:,3+(i-1)*3);
        z_s = s_e_n_z(:,3+(i-1)*3);
        ge_ = logical(keepI(1,i));
        gn_ = logical(keepI(2,i));
        gz_ = logical(keepI(3,i));
        %disp([ge_ gn_ gz_])
        
        subplot(3,3,1+3*(rowNum_-1)); hold on;
        if ge_
            plot(tdum,af*e,'linewidth',LineWidth); zoom on; grid on;
        else
            plot(tdum,af*e,'color',[0.5 0.5 0.5],'linewidth',LineWidth); zoom on; grid on;
        end
        title(strcat(stnms{i},', E'));
        subplot(3,3,2+3*(rowNum_-1)); hold on
        if gn_
            plot(tdum,af*n,'linewidth',LineWidth); zoom on; grid on;
        else
            plot(tdum,af*n,'color',[0.5 0.5 0.5],'linewidth',LineWidth); zoom on; grid on;
        end
        title(strcat(stnms{i},', N'));
        subplot(3,3,3+3*(rowNum_-1)); hold on;
        if gz_
            plot(tdum,zf*af*z,'linewidth',LineWidth); zoom on; grid on;
        else
            plot(tdum,zf*af*z,'color',[0.5 0.5 0.5],'linewidth',LineWidth); zoom on; grid on;
        end
        
        subplot(3,3,1+3*(rowNum_-1)); hold on;
        plot(tdum,af*e_s,'linewidth',LineWidth); zoom on; grid on;
        xlim([0 maxdur]);
        xlabel('Time [sec.]')
        subplot(3,3,2+3*(rowNum_-1)); hold on;
        plot(tdum,af*n_s,'linewidth',LineWidth); zoom on; grid on;
        xlim([0 maxdur]);
        xlabel('Time [sec.]')
        subplot(3,3,3+3*(rowNum_-1)); hold on;
        plot(tdum,zf*af*z_s,'linewidth',LineWidth); zoom on; grid on;
        xlim([0 maxdur]);
        xlabel('Time [sec.]');
        title(strcat(stnms{i},', Z'));
        text(maxdur-5,dist(i),stnms{i},'FontSize',14);
    end
end