function S = pickWaveforms(S,lfc,hfc,defaultP)
%
% pickWaveforms loops through traces and pick 1 phase and 1 polarity
%
% S = pickWaveforms(S,lfc,hfc,defaultP)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
if nargin < 2
    lfc = -inf;
end

if nargin < 3
    hfc = -inf;
end

if nargin < 4
    defaultP = [];
end

%%
sizeS = size(S);
refs = pull(S,'ref');
badI = isnat(refs);
goodEventsI = sum(~badI) == sizeS(1); % index to events (columns) with NO bad traces... if any bad traces exist, then the entire column (event) is bad
S(:,~goodEventsI) = []; % delete bad events (non-good events...)

refs = min(pull(S,'ref'))';
lEvents = sum(goodEventsI); %sizeS(2);
goodEventsI = true(lEvents,1);
button = 1;
n = 0;
while true
    if button == 2 || button == 98 % if middle click (2) or `b` button (98)
        if n == 1
            %we're assuming middle or b buttons were accidentally clicked because
            %if n == 1 we cant go to the previous trace (there is no previous trace)
            button = 1;
            n = 0;
            continue; % <-- go to 11
        end
        n = n - 1;
    else
        n = n + 1;
    end

    while true % plot waveforms and ready for ginput
        polarity = 0;
        close all;
        thisS = S(:,n);
        thisS = syncWaveforms(thisS,false,true);
        thisS = resampleWaveforms(thisS,200);

        fprintf("event number: %d, %s.%s.%s\n",n,thisS(1).knetwk,...
            thisS(1).kstnm,thisS(1).kcmpnm);

        [~,ax] = plotWaveforms(thisS,lfc,hfc,'k-',4,true);
        hold(ax,"on");
        ylims = ylim;
        xlims = xlim;

        if n > 1
            titleStr = sprintf("Current Trace: %d; Previous Trace: %d;" + ...
                "previous time: %f [sec.]; previous polarity: %d; keep?: %d",...
                n,n-1,S(1,n-1).user0,S(1,n-1).user2,goodEventsI(n-1));
            fprintf("%s\n",titleStr);
            ax(1).Title.String = titleStr;
        else
            ax(1).Title.String = sprintf("Current Trace: %d",n);
        end

        %% initialize
        if isempty(defaultP)
            x = mean(xlims);
        else
            if isfloat(defaultP)
                x = defaultP;
            end
        end

        h = plot(ax(1),repmat(x,2,1),ylims,'r--','linewidth',1);
        xlim(xlims);

        %% zoom to desired level so that you can make a pick, press any button to activate crosshairs
        pause(); % click any button to exit "pause" and activate crosshairs
        hold(ax,'on');
        xlims = xlim;
        ylims = ylim;

        %     default polarity is 0;
        %     1 = left-click; 2 = middle-click; 3 = right-click
        %     P = 112; S = 115;
        %     space-bar = 32; u = 117; d = 100; q == 113;
        %     o == 111; n = 110;

        %% pick
        while true % left click as many times as you want
            h.Visible = 'off';
            h = plot(repmat(x,2,1),ylims,'k--','linewidth',1);

            %%
            [x,y,button] = ginput(1);
            x = num2ruler(x(end),ax(1).XAxis); %converts x to a datetime
            y = y(end);
            button = button(end);

            %%
            fprintf("%f [sec.], %s, button: %d\n",seconds(x),refs(n)+x,button);
            h.Visible = 'off';
            if button ~= 1 %<--only when button OTHER than left click is clicked do we exit this loop, otherwise, keep recording ginput!
                break; % <-- go to 126
            end
        end

        %
        % ok, ive exited loop because ive clicked something other than left click (1)
        %

        %
        % if button clicked is "o" (111), then go to 51 (close and do over)
        %

        %%
        if button == 111 %|| button == 98 || button == 2
            continue;
        else
            % button other than "o" (111), middle click (2) or "b" (98) clicked
            % this means im happy with the P or S pick
            break; % <-- else, go to 136
        end
    end

    %%
    if button == 98 || button == 2
        % if button == 2 or 98 (middle click or b button), this means the user
        % made a mistake and would like to go back to the _previous_ trace.
        % this only work up until the very last trace.
        fprintf("\n");
        continue; % <-- go to 38
    end

    %% disp output (or mark for deletion)
    %fprintf("button: %d\n",button);
    if x < xlims(1) || x > xlims(2) || y < ylims(1) || y > ylims(2)
        %if last click is outside of window, ignore button and mark for deletion
        goodEventsI(n) = false;
        fprintf('marked trace %d for deletion\n',n);
    else
        fprintf('picked P phase at %s\n',refs(n) + x);
        goodEventsI(n) = true;

        if button == 117        % button is u for up (positive polarity)
            polarity = 1;
        elseif button == 100    % button is d for down (negative)
            polarity = -1;
        end
        S(1,n).user2 = polarity;

        if button == 115 % "S" clicked
            S(1,n).user1 = seconds(x);
        else
            S(1,n).user0 = seconds(x);
        end
        fprintf("\n");
    end

    %% exit while loop
    if n >= lEvents
        break;
    end
end

%%
% rI = find(rI);
% S = S(goodI);
% Sorig(rI(goodI)) = S;
% S = Sorig;