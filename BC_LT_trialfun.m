function [iv_inhb,iv_noInhb, iv_running] = BC_LT_trialfun(pos, iv_inhb, plot_flag)
%%
if nargin < 3
    plot_flag = 0
end
%%
  x_data=pos.data(3,:);
    %assign your running area
    max_tresh=75;
    min_tresh=20;
    x_running_index= find(x_data>=min_tresh & x_data<=max_tresh); %Find the data within your treshold
    x_dif=diff(x_running_index); 
    z=find(x_dif>1);%Gen idx of the gaps, calulated by the differential
    y=z+1; %Z has the values of the gaps beginnings
    jumps_end=x_running_index(z);
    jumps_end=[jumps_end,x_running_index(end)];
    jumps_start=x_running_index(y);
    jumps_start=[x_running_index(1), jumps_start];
    jump_start_time=pos.tvec(jumps_start);%Get the time of the start of the jumps from the tvec
    jump_end_time=pos.tvec(jumps_end);%Get the time of the end of the jumps from the tvec
    iv_running=iv(jump_start_time, jump_end_time); % Created the interval where the mouse is running from the ends and beginning
    % ----To do----
    %1.Remove those intervals where the mouse is to close to the
    %previous one
    %
    %2.Check for velocity
    %---
    
    % Substarct the iv where the mosue is inhibited from those where it is running
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
        surface([x;x],[y;y],[z;z],[col;col],... % This create a surface in the shape of a linear plot, this way it allows to mdify the color of the line accoridn to the speed of the mouse
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
        
        %Plot of position and inhibitions
        ax2=subplot(2,1,2);
        plot(pos.tvec,x_data, 'Color', BC_color_genertor('Oxford_blue'));
        h01=LTplotIvBars(iv_inhb,x_data,BC_color_genertor('Archt_green'),0.8);
        h02=LTplotIvBars(iv_noInhb,x_data,BC_color_genertor('Burnt_orange'),0.4);
        xlim([pos.tvec(1) pos.tvec(end)]);
        ylim([0 100]);
        % xlabel('Time (s)','Fontsize',14);
        % ylabel('Position in the maze (cm)');
        leg=legend([h01;h02],{'Light','No light'});
        legendFontSize = 12; % Adjust the font size as needed
        leg.Position=[0.890 0.3914 0.1047 0.0487];
        set(leg, 'FontSize', legendFontSize);
        box off;
        legend boxoff ;
        set(gca, 'TickDir', 'out');
        

        %% (*)Aesthetics (*)
        %Give common xlabel, ylabel, and title
        han=axes(fig,'visible','off');
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        han.YLabel.Position=[-0.0300 0.5000 0];
        ylabel(han,'Position in the x axis (cm)','FontSize', 16);
        han.XLabel.Position=[0.5000 -0.060 0];
        xlabel(han,'Time(s)','FontSize', 16);
        title(han,'Example of mouse displacement during linear track','FontSize', 20);
        fig = gcf;                   % Get current figure handle
        fig.Color = [1 1 1];
        fig.Color = [1 1 1];         % Set background color to white
        fig.Position = [100, 100, 1600, 700];  % [x, y, width, height]
        hold off;
        
    end

