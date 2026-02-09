function S = pickWaveformsMultiple(S,lfc,hfc,defaultP)
%
% pickWaveformsMultiple loops through traces and pick 1 phase and 1 polarity
%
% S = pickWaveformsMultiple(S,lfc,hfc,defaultP)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Sep 18, 2019

%%
if nargin < 2; lfc = -inf; end
if nargin < 3; hfc = -inf; end
if nargin < 4; defaultP = []; end

Sorig = S;
ref = pull(S(:,1),'ref');
rI = ~isnat(ref);
S = S(rI,:);
ref = ref(rI);

[lS,nc] = size(S);
goodI = true(lS,1);
button = 1;
n = 0;
x = NaT(nc,1);
polarity = zeros(1:nc);   % default polarity is 0;
while true
    if button == 2 || button == 98
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
        disp(n);
        close all;
        S(n,1:nc) = detrendWaveforms(syncWaveforms(S(n,1:nc)))';
        refTime = S(n,1).ref;
        ax = plotWaveforms(S(n,1:nc),lfc,hfc,'.-',4);
        
        xlims = NaT(nc,2);
        ylims = zeros(nc,2);
        if n > 1
            prevTime = ref(n-1) + S(n-1,currentLineAxesNumber).user0;
            if isnat(prevTime)
                suptitle({['Previous Trace: ',num2str(n-1),...
                    '; previous polarity: ',num2str(S(n-1).user1),...
                    '; keep?: ',num2str(goodI(n-1))];...
                    ' ';...
                    ['Current Trace: ',num2str(n)]});
            else
                suptitle({['Previous Trace: ',num2str(n-1),...
                    '; previous time: ',datestr(prevTime),...
                    '; previous polarity: ',num2str(S(n-1).user1),...
                    '; keep?: ',num2str(goodI(n-1))];...
                    ' ';...
                    ['Current Trace: ',num2str(n)]});
            end
        else
            suptitle(['Current Trace: ',num2str(n)]);
        end
        
        %% initialize
        if isempty(defaultP)
            xPlot(1:nc) = refTime + 0.5*(max(xlims) - min(xlims));
        else
            if isfloat(defaultP)
                xPlot(1:nc) = refTime + seconds(defaultP);
            end
        end
        
        for kk = 1:nc
            hold(ax(kk),'on');
            xlims(kk,:) = ax(kk).XLim;
            ylims(kk,:) = ax(kk).YLim;
            h(kk) = plot(ax(kk),repmat(xPlot(kk),2,1),ylims(kk,:),'r-','linewidth',2);
        end
        f = gcf;
        allAxes = findobj(f.Children,'Type','Axes');
        
        %% zoom to desired level so that you can make a pick, press any button to activate crosshairs
        pause(); %click any button to exit "pause" and activate crosshairs
        
        %     1 = left-click; 2 = middle-click; 3 = right-click
        %     P = 112; space-bar = 32; u = 117; d = 100; q == 113;
        %     o == 111; n = 110;
        
        %% pick
        while true % left click as many times as you want
            for kk = 1:nc
                delete(h(kk));
                ax(kk).ColorOrderIndex = 3;
                if ~isnat(x(kk))
                    h(kk) = plot(ax(kk),repmat(x(kk),2,1),ylims(kk,:),'-','linewidth',2);
                else
                    h(kk) = plot(ax(kk),repmat(xPlot(kk),2,1),ylims(kk,:),'-','linewidth',2);
                end
            end
            
            %%
            [xTmp,yTmp,button] = ginput(1);
            currentAxis = gca;
            xTmp = num2ruler(xTmp,currentAxis.XAxis);         %converts x to a datetime
            numClicked = find(currentAxis == allAxes);
            allAxes(numClicked).ColorOrderIndex = 3;
            hclick = plot(allAxes(numClicked),repmat(xTmp,2,1),allAxes(numClicked).YLim,'-','linewidth',2);
            currentLineAxesNumber = numClicked;
            if currentLineAxesNumber > 2
                currentLineAxesNumber = currentLineAxesNumber - 1;
            end
            
            %%
            x(currentLineAxesNumber) = xTmp;
            if button == 117        % button is u for up (positive polarity)
                polarity(currentLineAxesNumber) = 1;
            elseif button == 100    % button is d for down (negative)
                polarity(currentLineAxesNumber) = -1;
            end
            
            if button == 110 || button == 98 || button == 2 || button == 111
                break; % <-- go to 83
            end
            hclick.Visible = 'off';
        end
        
        if button ~= 111 % if button == 111, go to 25
            % o = 111;
            % for any button other than 'o,' we save x, y, and button and
            % exit loop
            % if its 'o' we replot _current_ trace
            break; % <-- else, go to 93
        end
    end
    
    %%
    if button == 98 || button == 2
        % if button == 2 or 98 (middle click or b button), this means the user
        % made a mistake and would like to go back to the _previous_ trace.
        % this only work up until the very last trace.
        continue; % <-- go to 11
    end
    
    %% disp output (or mark for deletion)
    disp(['button: ',num2str(button)])
    xlims = xlims(1,:);
    ylims = ylims(1,:);
    if xTmp < xlims(1) || xTmp > xlims(2) || yTmp < ylims(1) || yTmp > ylims(2)
        %if last click is outside of window, ignore button and mark for deletion
        goodI(n) = false;
        disp(['marked trace ',num2str(n),'for deletion']);
    else
        for kk = 1:nc
            if ~isnat(x(kk))
                fprintf('picked P phase at %s:\n',datestr(x(kk)));
                goodI(n) = true;
                S(n,kk).user0 = x(kk)-refTime;
                S(n,kk).user1 = polarity(kk);
            end
        end
    end
    
    %% exit while loop
    if n >= lS
        break;
    end
end

%%
rI = find(rI);
S = S(goodI,:);
Sorig(rI(goodI),:) = S;
S = Sorig;
