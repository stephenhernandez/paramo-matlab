% Code for picking arrivals on traces

%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:   Morgan Plain %
% DATE:     18/11/2015   %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all

% load data
%load('Data_random/morgan_test_events.mat')
events = indiv_events;

% note: input data is in the form of a matrix (called 'events') where
%       each column is the data for one trace

% settings
fullFlag = true;     %(0=normal, 1=fullscreen)
zoomFlag = false;           %(0=normal, 1=new zoom window)


%% format inputs
nTraces = size(events,2); %how many traces are there

% traces
for i = 1:nTraces %loop over traces
    traces{i}(:,1) = events(:,i); %amplitude
end

% picks
for i = 1:nTraces
    % picks(i) = dataset_picks(i);  % use this for variable picks
    % single number cell - time of pick for the corresponding trace
    picks(i) = 400; %constant at 400
end
% note: ASSUMES ALL TRACES SAME LENGTH


%% instructions
fprintf('\nTRACE PICKER\n')
fprintf(' - original pick plotted in red\n')
fprintf(' - if zoom is on, left click to choose zoom centre-point\n')
fprintf(' - left click to make a new pick \n')
fprintf(' - new pick plots in green\n')
fprintf(' - you can pick as many times as you like\n')
fprintf(' - confirm the last pick by right clicking\n')
fprintf(' - moves to next trace\n\n')


%% plot

% setting variables
delcount=0; % for_deletion counter

% loop over all traces
for trace_num = 1:nTraces
    
    % announce trace number
    fprintf('\nCurrently picking trace #%d/%d \n',trace_num,nTraces)
    
    % plot trace
    hfig = figure(100);
    if fullFlag == 1 %make figures fullscreen
        set(hfig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    plot(traces{trace_num}(:,1))
    hold on
    % make line vector and plot original pick
    ylims = ylim;
    xlims = xlim;
    picky(:,1) = linspace(ylims(1),ylims(2),10)';
    picky(:,2) = (ones(length(picky),1)*picks(trace_num))';
    plot(picky(:,2),picky(:,1),':','linewidth',3); %plot initial pick right here with red line
    
    %title(['Event: ',num2str(trace_num),'/',num2str(trace_num1),'; ',datestr(twin(trace_num))])
    % plot zoom window
    if zoomFlag %create zoom window
        [pick_time,y] = ginput(1); % user clicks zoom centre point
        pick_time = round(pick_time);
        l1 = pick_time-100;
        if l1 < 0 % stop an outside-of-axis error
            l1 = 0;
        end
        
        xzoomlim(1,1:2) = [l1 pick_time+100]; % set zoom limits
        hfig = figure(101);
        plot(xzoomlim(1):xzoomlim(2),traces{trace_num}(xzoomlim(1):xzoomlim(2)),'b','linewidth',3)
        title('Zoomed trace')
        if picks(trace_num)>xzoomlim(1) && picks(trace_num)<xzoomlim(2) % only plot old pick if within zoom axis
            hold on
            picky(:,1) = linspace(ylims(1),ylims(2),10)';
            picky(:,2) = (ones(length(picky),1)*picks(trace_num))';
            plot(picky(:,2),picky(:,1),'r:','linewidth',3);
        end
        
        % change xlim variable to be new zoomed limits
        % (so that marking for deletion works)
        xlims = xzoomlim;
    end %if zoom
    
    % resetting for new trace
    button = 1; %1=left click
    clear pick_time %delete the variable
    
    % make a new pick
    while button == 1
        a = exist('pick_time','var'); %is this the first iteration (0=yes)
        [pick_time,y,button] = ginput(1);
        
        % if click is outside figure, mark for deletion
        if pick_time<xlims(1) || pick_time>xlims(2) || y<ylims(1) || y>ylims(2)
            delcount=delcount+1;
            for_deletion(delcount,1) = trace_num;
            fprintf('Marked for deletion. \n')
            break %leave while loop and go to next trace
        end
        
        if button == 1
            fprintf('Picked time: %f \n', pick_time);
            final_picks(trace_num,1) = pick_time;
            final_picks(trace_num,2) = abs(pick_time-picks(trace_num));
        elseif button == 3 %right click
            fprintf('Confirming pick..\n')
        end
        
        % if right clicked, use last pick (unless it's the first iteration)
        if button == 3
            if a == 1 %if not first iteration
                fprintf('Confirmed pick at time: %f for trace #%d \n',final_picks(trace_num,1),trace_num)
                pause(1);
                continue %finish while loop -> next trace
            end
            fprintf('  Error: You did not chose a pick before confirming\n')
            fprintf('  You can either confirm the last pick (Rclick) or pick again\n\n')
        end
        
        %plot new pick
        picknew(:,1) = linspace(ylims(1),ylims(2),10)';
        picknew(:,2) = ones(length(picknew),1)*pick_time;
        hold on
        plot(picknew(:,2),picknew(:,1),'g:','linewidth',3);
        
    end %while
    
    % close trace plot
    close(100)
    if zoomFlag
        close(101)
    end
    
end %for trace_num

fprintf('\n\nAll traces picked.\n')
try
    fprintf('  %d traces marked for deletion.\n',length(for_deletion))
end
