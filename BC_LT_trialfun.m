function [iv_inhb,iv_noInhb, iv_running] = BC_LT_trialfun(pos, iv_inhb, plot_flag)
%%

Archt_green = [0.530 0.820 0.645];
Oxford_blue = [0.039 0.137 0.259];
Powder_blue = [0.682 0.772 0.921];
Burnt_orange= [0.758 0.348 0.249];
Red_crayola= [0.937 0.176 0.337];
Web_orange= [1.00 0.678 0.020];

Archt_green_alpha = [0.530 0.820 0.645 0.3];
Oxford_blue_alpha= [0.039, 0.137, 0.259, 0.300];
Powder_blue_alpha= [0.682 0.772 0.921 0.3];
Burnt_orange_alpha= [0.933 0.423 0.302 0.3];
Web_orange_alpha= [1.00 0.678 0.020 0.3];

if nargin < 3
    plot_flag = 1
end

%%
  x_data=pos.data(3,:);
    %assign your running area
    max_tresh=75;
    min_tresh=20;
    x_running_index= find(x_data>=min_tresh & x_data<=max_tresh);
    x_dif=diff(x_running_index);
    z=find(x_dif>1);
    y=z+1;
    jumps_end=x_running_index(z);
    jumps_end=[jumps_end,x_running_index(end)];
    jumps_start=x_running_index(y);
    jumps_start=[x_running_index(1), jumps_start];
    % f_jumps=[jumps_start,jumps_end];
    % sort(f_jumps);
    % %plot thhe data
    % plot(x_data)
    % hold on
    % scatter(jumps_end,x_data(jumps_end),'filled')
    % scatter(jumps_start,x_data(jumps_start),'filled')
    %Convert this indices to time epochs
    jump_start_time=pos.tvec(jumps_start);
    jump_end_time=pos.tvec(jumps_end);
    iv_running=iv(jump_start_time, jump_end_time); % This is the interval where the mouse is running
    % ----To do----
    %1.Remove those intervals where the mouse is to close to the
    %previous one
    %
    %2.Check for velocity
    %---
    
    % Substarct the iv where the mosue is inhibited from those where it is
    % running
    cfg01=[];
    cfg01.verbose=0;
    iv_noInhb=DifferenceIV(cfg01,iv_running,iv_inhb);
    
    
    if plot_flag
        % Plot the x_pos by time and plot the int to collaborate that you got them right
        fig=figure(1919);
        clf;
        x = pos.tvec;
        y = x_data;
        z = zeros(size(x));
        speed_data=pos.data(5,:);
        col = speed_data;  % This is the color, it varies with x in this case.
        ax1=subplot(2,1,1);
        surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2);
        c=colorbar;
        c.Label.String='Velocity (cm/s)';
        c.Label.FontSize=12;
        c.Position= [0.910 0.5900 0.0100 0.330];
        xlim([pos.tvec(1) pos.tvec(end)]);
        set(gca, 'TickDir', 'out');
        
        ax1=subplot(2,1,1);
        
        %xlabel('Time (s)');
        %ylabel('Position in the x axis of the maze (cm)');
        %h03=LTplotIvBars(iv_running,x_data,Oxford_blue,0.3)
        %corrected_time=csc.tvec-min(csc.tvec); % creates n array of time corrected for the time that neuralynx started recording
        ax2=subplot(2,1,2);
        plot(pos.tvec,x_data, 'Color', Oxford_blue);
        h01=LTplotIvBars(iv_inhb,x_data,Archt_green,0.8);
        h02=LTplotIvBars(iv_noInhb,x_data,Burnt_orange,0.4);
        xlim([pos.tvec(1) pos.tvec(end)]);
        ylim([0 100]);
        % xlabel('Time (s)','Fontsize',14);
        % ylabel('Position in the maze (cm)');
        leg=legend([h01;h02],{'Laser 520nm','Running no inhibition'});
        legendFontSize = 12; % Adjust the font size as needed
        leg.Position=[0.850 0.3914 0.1047 0.0487];
        set(leg, 'FontSize', legendFontSize);
        box off;
        legend boxoff ;
        set(gca, 'TickDir', 'out');
        
        %Give common xlabel, ylabel, and title
        han=axes(fig,'visible','off');
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        han.YLabel.Position=[-0.0300 0.5000 0];
        ylabel(han,'Position in the x axis of the maze (cm)','FontSize', 16);
        han.XLabel.Position=[0.5000 -0.060 0];
        xlabel(han,'Time(s)','FontSize', 16);
        title(han,'Example of displacement of mouse in the x axis during inhibition session','FontSize', 20);
        fig = gcf;                   % Get current figure handle
        fig.Color = [1 1 1];
        hold off;
        
    end

