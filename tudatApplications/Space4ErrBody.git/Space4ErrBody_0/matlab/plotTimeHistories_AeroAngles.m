function [  ] = plotTimeHistories_AeroAngles( compilation, mainpath )


%% Time History: Angle of Attack - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 6000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack  through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-10 50])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', -10:5:50);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack);
        end
        
        %plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/AoA_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
         close(fig_num);
    end
end


%% Time History: Evaluated Angle of Attack - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 6000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Angle of Attack  through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-10 50])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Evaluated Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', -10:5:50);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_angle_of_attack);
        end
        
        %plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/Evaluated_AoA_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
         close(fig_num);
    end
end

%% Time History: Bank Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 6000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle  through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-90 90])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
        end
        
        %plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/bank_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
         close(fig_num);
    end
end

%% Time History: Evaluated Bank Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 6000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle  through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Evaluated Bank Angle (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_bank_angle);
        end
        
        %plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/evaluated_bank_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
         close(fig_num);
    end
end



%% Time History: Flight-Path Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 654000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Flight-Path Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-90 90])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Flight-Path Angle (deg)') % y-axis label
        set(gca,'YTick', -90:30:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.flight_path_angle);
        end
        
        %plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/flight_path_angle_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



end