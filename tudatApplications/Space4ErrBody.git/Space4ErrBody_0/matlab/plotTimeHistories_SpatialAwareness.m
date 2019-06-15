function [  ] = plotTimeHistories_SpatialAwareness( compilation )

%% Time History: Cumulative Cartesian Distance Travelled - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 833000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat(' Cumulative Cartesian Distance Travelled through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 8000])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel(' Cumulative Cartesian Distance Travelled (km)') % y-axis label
        set(gca,'YTick', 0:500:8000);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.cumulative_distance_travelled/1e3);
        end
        %plot([0 4000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/cumCartDisTravelled_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end



%% Time History: Height - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 823000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        
        ylim([10 130])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 10:10:130);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],(10)*[1 1],'k','LineWidth',2)
        
        
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,'k','LineWidth',2);
                            xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([20 130])
                set(gca,'YTick', 20:10:130);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
            end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end

%% Time History: Distance Traveled - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Distance Traveled through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 180])
        max_tof = max([compilation(p).evolutions.max_tof]);
        %xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Distance Traveled (deg)') % y-axis label
        %set(gca,'YTick', 0:30:180);
        %set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.distance_traveled);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/distance_traveled_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end

%% Time History: Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Distance To Go through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 60])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Distance To Go (deg)') % y-axis label
        set(gca,'YTick', 0:10:60);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],(.75)*[1 1],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,'k','LineWidth',2);
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                %ylim([0 1])
                %set(gca,'YTick', 0:.2:2);
                plot([0 1400],(.75)*[1 1],'k','LineWidth',2)
                
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go);
            end
        end
        
        
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/distance_to_go_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end

%% Time History: Distance Covered Ratio - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Distance Covered Ratio - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Distance Covered Ratio (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],(.4)*[1 1],'k','LineWidth',2)
        plot([0 max_tof],(.1)*[1 1],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if compilation(p).evolutions(k).trajectories(ii).individual.time_vector(end) > 1660
                distanceToCover = compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(1);
                distanceCoveredRatio = compilation(p).evolutions(k).trajectories(ii).individual.distance_traveled/distanceToCover;
                h = stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    distanceCoveredRatio);
                   set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))

            terminationDistanceRatio = compilation(p).evolutions(k).printedPopulationDV.terminationDistanceRatio(ii);
                   
                I = find( distanceCoveredRatio > terminationDistanceRatio );
   scatter(compilation(p).evolutions(k).trajectories(ii).individual.time_vector(I(1)),terminationDistanceRatio);
    'e';
                   
            end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/distanceCoveredRatio_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end


%% Time History: Heading Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 783000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 360])
        xlim([0 1400])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading Angle (deg)') % y-axis label
        set(gca,'YTick', 0:30:360);
        set(gca,'XTick', 0:200:1400);
        hold on
        grid on
        %  for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.heading_angle);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heading_angle_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Time History: Heading to Target - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 723000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading to Target through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 360])
        xlim([0 4000])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading to Target (deg)') % y-axis label
        set(gca,'YTick', 0:30:360);
        set(gca,'XTick', 0:200:4000);
        hold on
        
        %  for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.heading_to_target);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heading_to_target_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Time History: Heading Error - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 734000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Error through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-360 360])
        xlim([0 4000])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading Error (deg)') % y-axis label
        set(gca,'YTick', -360:30:360);
        set(gca,'XTick', 0:200:4000);
        hold on
        grid on
        %  for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.heading_error,'k','LineWidth',2);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heading_error_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
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
% for ii = 1:numel(compilation(p).evolutions(k).trajectories)
%             plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
%                 compilation(p).evolutions(k).trajectories(ii).individual.altitude/1e3);
%         end
%
%         plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
%         hold off
%         saveas(...
%             figure(fig_num),...
%             strcat(...
%     compilation(p).mainpath,...
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
% for ii = 1:numel(compilation(p).evolutions(k).trajectories)
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