function [  ] = plotTimeHistories_SpatialAwareness( compilation, mainpath )

%% Time History: Height - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 823000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 150])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:10:150);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/height_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Time History: Distance Traveled - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Distance Traveled through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 180])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Distance Traveled (deg)') % y-axis label
        set(gca,'YTick', 0:30:180);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.distance_traveled);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/distance_traveled_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Time History: Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Distance To Go through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 180])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Distance To Go (deg)') % y-axis label
        set(gca,'YTick', 0:30:180);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/distance_to_go_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Time History: Heading to Target - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading to Target through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-180 180])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading to Target (deg)') % y-axis label
        set(gca,'YTick', -180:30:180);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.heading_to_target);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/heading_to_target_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Time History: Heading Error - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Error through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-60 60])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading Error (deg)') % y-axis label
        set(gca,'YTick', -60:15:60);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.heading_error);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/heading_error_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



% %% Time History: Radial Distance - per Evolution
% for p = 1:numel(compilation)
%     
%     for k = 1:numel(compilation(p).evolutions)
%         fig_num = p*100 + 6000 + k*1;
%         figure(fig_num)
%         set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
%         set (gca,'Fontsize',15)
%         title(strcat('R_{norm} through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
%         ylim([6.35e3 6.6e3])
%         max_tof = max([compilation(p).evolutions.max_tof]);
%         xlim([0 max_tof])
%         xlabel('Propagation Time (s)') % x-axis label
%         ylabel('Radial Distance (km)') % y-axis label
%         set(gca,'YTick', 6.350e3:0.05e3:6.5e3);
%         set(gca,'XTick', 0:200:max_tof);
%         hold on
%         
%         for ii = 1:numel(compilation(p).evolutions(k).trajectories)
%             plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
%                 compilation(p).evolutions(k).trajectories(ii).individual.altitude/1e3);
%         end
%         
%         plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
%         hold off
%         saveas(...
%             figure(fig_num),...
%             strcat(...
%             mainpath,...
%             '/figures/altitude_T_Evolution_',...
%             num2str(k - 1),...
%             '_Set',...
%             convertCharsToStrings(compilation(p).set),...
%             '.png'),...
%             'png');
%         close(fig_num);
%     end
% end
%
%
% %% Time History: Radial Distance - per Set
% for p = 1:numel(compilation)
%
%     fig_num = p + 7000;
%     figure(fig_num)
%     set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
%     set (gca,'Fontsize',15)
%     title(strcat('R_{norm} through Time - per Set ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
%     ylim([6.35e3 6.6e3])
%     max_tof = max([compilation(p).evolutions.max_tof]);
%     xlim([0 max_tof])
%     xlabel('Propagation Time (s)') % x-axis label
%     ylabel('Radial Distance (km)') % y-axis label
%     set(gca,'YTick', 6.350e3:0.05e3:6.6e3);
%     set(gca,'XTick', 0:200:max_tof);
%     hold on
%
%     for k = 1:numel(compilation(p).evolutions)
%
%         for ii = 1:numel(compilation(p).evolutions(k).trajectories)
%             plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
%                 compilation(p).evolutions(k).trajectories(ii).individual.R_R_norm/1e3);
%         end
%     end
%
%     plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
%     hold off
%     saveas(figure(fig_num),strcat(mainpath,'/figures/R_norm_v_T_Set',...
%         convertCharsToStrings(compilation(p).set),...
%         '.png'),...
%         'png');
%     close(fig_num);
%
% end
%

end