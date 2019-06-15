function [  ] = plotTimeHistories_EngineOps( compilation )

%% Time History: Evaluated Throttle Setting - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 722000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Throttle Setting through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Evaluated Throttle Setting (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
for ii = 1:numel(compilation(p).evolutions(k).trajectories)   
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_throttle_setting);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluatedThrottleSetting_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end


%% Time History: Commanded Throttle Setting - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Commanded Throttle Setting through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Throttle Setting (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
for ii = 1:numel(compilation(p).evolutions(k).trajectories)   
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.commanded_throttle_setting);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/commandedThrottleSetting_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end


%% Time History: Engine Status - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Engine Status through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 2])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Engine Status (-)') % y-axis label
        set(gca,'YTick', 0:1:2);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
for ii = 1:numel(compilation(p).evolutions(k).trajectories)   
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.engine_status);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/engine_status_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %        close(fig_num);
    end
end




%% Time History: Evaluated Thrust Elevation Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Thrust Elevation Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Evaluated Thrust Elevation Angle (deg)') % y-axis label
        set(gca,'YTick', -30:5:30);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
            % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_thrust_elevation_angle);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluated_thrust_elevation_angle_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end

%% Time History: Commanded Thrust Elevation Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Commanded Thrust Elevation Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Thrust Elevation Angle (deg)') % y-axis label
        set(gca,'YTick', -30:5:30);
        set(gca,'XTick', 0:200:max_tof);
        hold on
                grid on
        
                        plot([0 max_tof],(0)*[1 1],'k','LineWidth',2)

                for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.commanded_thrust_elevation_angle);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/commanded_thrust_elevation_angle_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end



%% Time History: Evaluated Thrust Azimuth Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Thrust Azimuth Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        max_tof = 2000;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Evaluated Thrust Azimuth Angle (deg)') % y-axis label
        set(gca,'YTick', -30:5:30);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
            % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_thrust_azimuth_angle);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluated_thrust_azimuth_angle_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %   close(fig_num);
    end
end


%% Time History: Commanded Thrust Azimuth Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Commanded Thrust Azimuth Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        max_tof = 2000;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Thrust Azimuth Angle (deg)') % y-axis label
        set(gca,'YTick', -30:5:30);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
            % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.commanded_thrust_azimuth_angle);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/commanded_thrust_azimuth_angle_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %   close(fig_num);
    end
end



end