function [  ] = plotCompoundRelations( compilation )



%% Time History: Height & Commanded Throttle Setting & Total g-Load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 823000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        
        ylim([0 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            
            h = stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/max(compilation(p).evolutions(k).trajectories(ii).individual.height))
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
            
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.commanded_throttle_setting);
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag/max(compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag));
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack/max(compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack));
            
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








%% Height vs. Freestream Temperature - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Freestream Temperature - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 130])
        xlim([150 400])
        xlabel('Temperature (K)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:10:130);
        set(gca,'XTick', 150:25:400);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.freestreamTemperature,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,'k','LineWidth',2);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/freestreamTemperature_v_height_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Angle of Attack vs. Mach Number - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack vs. Mach Number - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 30])
        xlabel('Mach Number (-)') % x-axis label
        ylabel('Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.mach,...
                    compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack,'k','LineWidth',2);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.mach,...
                    compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack);
            end
        end
        plot(compilation(p).evolutions(k).trajectories(ii).individual.alphaMachBounds_Mach,compilation(p).evolutions(k).trajectories(ii).individual.alphaMachBounds_UB,'k','LineWidth',2);
        plot(compilation(p).evolutions(k).trajectories(ii).individual.alphaMachBounds_Mach,compilation(p).evolutions(k).trajectories(ii).individual.alphaMachBounds_LB,'k','LineWidth',2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleOfAttack_v_mach_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Skip Suppression Limit vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %  for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    %for k = 30:numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 101000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Skip Suppression Limit vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(rad2deg([-pi/2 pi/2]))
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Skip Suppression Limit (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        plot([0 60],[0 0],'k','LineWidth',2)
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.skip_suppression_limit);
            %  stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %     compilation(p).evolutions(k).trajectories(ii).individual.evaluated_bank_angle);
            % stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
            %
            %                 stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %                     compilation(p).evolutions(k).trajectories(ii).individual.commanded_bank_angle);
            
            
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/skip_suppression_limit_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end

%% Skip Suppression Limit vs. Normalized Specific Energy - per Evolution
for p = 1:numel(compilation)
    
    %  for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    %for k = 30:numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 101000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Skip Suppression Limit vs. Normalized Specific Energy - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(rad2deg([-pi/2 pi/2]))
        xlim([0 1])
        xlabel('Normalized Specific Energy (-)') % x-axis label
        ylabel('Skip Suppression Limit (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:.1:1);
        hold on
        grid on
        plot([0 60],[0 0],'k','LineWidth',2)
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                compilation(p).evolutions(k).trajectories(ii).individual.skip_suppression_limit);
            %  stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %     compilation(p).evolutions(k).trajectories(ii).individual.evaluated_bank_angle);
            % stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
            %
            %                 stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %                     compilation(p).evolutions(k).trajectories(ii).individual.commanded_bank_angle);
            
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/skip_suppression_limit_v_normalizedSpecificEnergy_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end


%% Total Passenger Z-Component G-load vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 20
        fig_num = p*100 + k*1 + 102000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Passenger Z-Component G-load vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-5 5])
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Total Passenger Z-Component G-load (g)') % y-axis label
        set(gca,'YTick', -5:1:5);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            %   if sum(abs(compilation(p).evolutions(k).trajectories(ii).individual.heading_error) > 30) == 0
            %       if compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 30
            %          if sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) > 2
            %             if compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy(end) < 0.3
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.passenger_fixed_total_g_load_z);
            %     end
            % end
            %    end
            %end
        end
        
        plot([0 60],(0)*[1 1],'k','LineWidth',2)
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_passenger_g_load_z_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Normalized Specific Energy vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %   for k = 30
        fig_num = p*100 + k*1 + 103000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Normalized Specific Energy vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1.5])
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Normalized Specific Energy (-)') % y-axis label
        set(gca,'YTick', 0:.1:1.5);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        plot([0 60],(1)*[1 1],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/normalizedSpecificEnergy_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Thrust Magnitude vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 104000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Magnitude vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 1100])
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Thrust Magnitude (kN)') % y-axis label
        %set(gca,'YTick', 0:100:1100);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.thrustMagnitude/1e3);
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustMagnitude_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Thrust Magnitude vs. Normalized Specific Energy - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 104000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Magnitude vs. Normalized Specific Energy - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 1100])
        xlim([0 1])
        xlabel('Normalized Specific Energy ()') % x-axis label
        ylabel('Thrust Magnitude (kN)') % y-axis label
        %set(gca,'YTick', 0:100:1100);
        set(gca,'XTick', 0:.1:1);
        hold on
        grid on
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                compilation(p).evolutions(k).trajectories(ii).individual.thrustMagnitude/1e3);
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustMagnitude_v_normalizedSpecificEnergy_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Thrust Acceleration Components vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 105000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Acceleration Components vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-40 40])
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Thrust Acceleration Components (m/s^2)') % y-axis label
        set(gca,'YTick', -40:10:40);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_x,'k');
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_y,'r');
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_z,'b');
            
        end
        
        legend('x-dir','y-dir','z-dir')
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAccelerationsComponents_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Thrust Acceleration Components vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 106000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Acceleration Magnitude vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 40])
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Thrust Acceleration Magnitude (m/s^2)') % y-axis label
        set(gca,'YTick', 0:5:40);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_M);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAccelerationMagnitude_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Height vs. Density - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 107000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Density - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 100])
        xlim([0 .05])
        xlabel('Density (kg/m^3)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:10:100);
        set(gca,'XTick', 0:.01:.05);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.localDensity,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,'k','LineWidth',2);
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_localDensity_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Density vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 108000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Density vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 800])
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Density (kg/m^3)') % y-axis label
        % set(gca,'YTick', 0:100:800);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            %   if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 20 )
            %     if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.localDensity);
            % end
            %     end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/density_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Height vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 109000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 150])
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Height (km)') % y-axis label
        % set(gca,'YTick', 0:25:150);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Mass vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 110000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Mass vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 150])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Mass (10^3 kg)') % y-axis label
        %set(gca,'YTick', 0:15:150);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.mass/1e3);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/mass_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Mass vs. Normalized Specific Energy - per Evolution
for p = 1:numel(compilation)
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 110000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Mass vs. Normalized Specific Energy - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 150])
        xlim([0 1])
        xlabel('Normalized Specific Energy (-)') % x-axis label
        ylabel('Mass (10^3 kg)') % y-axis label
        %set(gca,'YTick', 0:15:150);
        set(gca,'XTick', 0:.1:1);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                compilation(p).evolutions(k).trajectories(ii).individual.mass/1e3);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/mass_v_normalizedSpecificEnergy_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Mass Rate vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 111000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Mass Rate vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 250])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Mass Rate (kg/s)') % y-axis label
        set(gca,'YTick', 0:25:250);
        set(gca,'XTick', 0:10:60);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                abs(compilation(p).evolutions(k).trajectories(ii).individual.mass_rate) );
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/massRate_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Mass vs. Height - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 112000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Mass vs. Height - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 150])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        % xlim([0 8000])
        xlabel('Height (km)') % x-axis label
        ylabel('Mass (kg)') % y-axis label
        % set(gca,'YTick', 0:10:150);
        % set(gca,'XTick', 0:500:8000);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
                compilation(p).evolutions(k).trajectories(ii).individual.mass);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/mass_v_height_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Height vs. Evaluated Throttle - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 113000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Throttle Setting vs. Height  - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 150])
        % xlim([0 8000])
        ylabel('Height (km)') % x-axis label
        xlabel('Evaluated Throttle Setting (-)') % y-axis label
        % set(gca,'YTick', 0:10:150);
        % set(gca,'XTick', 0:500:8000);
        hold on
        grid on
        %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.evaluated_throttle_setting,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluatedThrottle_v_height_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Evaluated Throttle Setting vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 114000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Throttle Setting vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([180 360])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Evaluated Throttle Setting (-)') % y-axis label
        %set(gca,'YTick', 180:30:360);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 30 )
                %    if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.evaluated_throttle_setting);
                %     end
            end
            
            
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluatedThrottleSetting_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Commanded Throttle Setting vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 115000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Commanded Throttle Setting vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([180 360])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Commanded Throttle Setting (-)') % y-axis label
        %set(gca,'YTick', 180:30:360);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.commanded_throttle_setting);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/commandedThrottleSetting_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Commanded Throttle Setting vs. Normalized Specific Energy - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 115000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Commanded Throttle Setting vs. Normalized Specific Energy - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([180 360])
        xlim([0 1])
        xlabel('Normalized Specific Energy (-)') % x-axis label
        ylabel('Commanded Throttle Setting (-)') % y-axis label
        %set(gca,'YTick', 180:30:360);
        set(gca,'XTick', 0:.1:1);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 30 )
                %    if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).trajectories(ii).individual.commanded_throttle_setting);
                %     end
            end
            
            
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/commandedThrottleSetting_v_normalizedSpecificEnergy_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Height vs. Normalized Specific Energy - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 116000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Normalized Specific Energy - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 130])
        xlim([0 1])
        xlabel('Normalized Specific Energy (-)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:10:130);
        set(gca,'XTick', 0:.1:1);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        %plot(V,compilation(1).evolutions(1).trajectories(1).individual.height/1e3)
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_normalizedSpecificEnergy_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end
%% Airspeed vs. Normalized Specific Energy- per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 116000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Airspeed va. Normalized Specific Energy - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 8000])
        xlim([0 1])
        ylabel('Airspeed (m/s)') % x-axis label
        xlabel('Normalized Specific Energy (-)') % y-axis label
        set(gca,'YTick', 0:500:8000);
        set(gca,'XTick', 0:.1:1);
        hold on
        grid on
        
        %  plot([0 1],(25)*[1 1],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                compilation(p).evolutions(k).trajectories(ii).individual.airspeed);
        end
        
        %plot(V,compilation(1).evolutions(1).trajectories(1).individual.height/1e3)
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/airspeed_v_normalizedSpecificEnergy_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Height vs. Airspeed - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 116000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Airspeed - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 130])
        xlim([0 8000])
        xlabel('Airspeed (m/s)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:10:130);
        set(gca,'XTick', 0:500:8000);
        hold on
        grid on
        
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.airspeed,...
                    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,'k','LineWidth',2);
               %             xlim([0 1400])
               % set(gca,'XTick', 0:200:1400);
                %ylim([0 1])
                %set(gca,'YTick', 0:.2:2);
              %  plot([0 1400],(.75)*[1 1],'k','LineWidth',2)
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.airspeed,...
                    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
            end
        end
        %plot(V,compilation(1).evolutions(1).trajectories(1).individual.height/1e3)
        
        plot([0 8000],(10)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_V_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Heading Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1
        fig_num = p*100 + k*1 + 117000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([180 360])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Heading Angle (deg)') % y-axis label
        set(gca,'YTick', 180:30:360);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        
        %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            %if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 20 )
            %   if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.heading_angle);
            %  end
            % end
            
            
        end
        
        %                 for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        %             %for ii = 1:numel(compilation(p).evolutions(k).trajectories)
        %            plot(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
        %                compilation(p).evolutions(k).trajectories(ii).individual.heading_to_target);
        %
        %         end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/headingAngle_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Heading Error vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 118000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Error vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(30*[-1 1])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Heading Error (deg)') % y-axis label
        set(gca,'YTick', -30:10:30);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        hold on
        grid on
        
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.heading_error,'k','LineWidth',2);
              
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.heading_error);
            end
        end
        
        
        
        
        plot(compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_distance,...
            compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_LB,'k','LineWidth',2);
        plot(compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_distance,...
            compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_UP,'k','LineWidth',2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/headingError_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end




%% Angle of Attack vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 119000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleOfAttack_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Angle of Attack vs. Distance Covered Ratio - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 119000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack vs. Distance Covered Ratio - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 1])
        xlabel('Distance Covered Ratio (-)') % x-axis label
        ylabel('Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        set(gca,'XTick', 0:.1:1);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            
            distanceToCover = compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(1);
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_traveled/distanceToCover,...
                compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleOfAttack_v_distanceCoveredRatio_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Angle of Attack vs. Normalized Specific Energy - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 119000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack vs. Normalized Specific Energy - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 1])
        xlabel('Normalized Specific Energy (-)') % x-axis label
        ylabel('Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        set(gca,'XTick', 0:.1:1);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
        
            h = stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack);
            %{
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
           
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_angle_of_attack);
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
            
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_angle_of_attack_Ascent);
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
            
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_angle_of_attack_Ascent);
                              set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
   
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_angle_of_attack_Descent);
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
            
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_angle_of_attack_Descent);
                     
            %}
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleOfAttack_v_normalizedSPecificEnergy_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end




%% Bank Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 120000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(90*[-1 1])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        hold on
        grid on
        
        plot([0 60],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngle_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Bank Angle vs. Normalized Specific Energy - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 120000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle vs. Normalized Specific Energy - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(90*[-1 1])
        xlim([0 1])
        xlabel('Normalized Specific Energy (-)') % x-axis label
        ylabel('Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:.1:1);
        % set(gca, 'XDir','reverse')
        hold on
        grid on
        
        plot([0 1],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngle_v_normalizedSpecificEnergy_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Evaluated Bank Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 121000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Bank Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(90*[-1 1])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Evaluated Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        plot([0 60],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_bank_angle);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluatedBankAngle_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Evaluated Bank Angle vs. Normalized Specific Energy - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 121000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Bank Angle vs. Normalized Specific Energy - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(90*[-1 1])
        xlim([0 1])
        xlabel('Normalized Specific Energy (-)') % x-axis label
        ylabel('Evaluated Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:.1:1);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        plot([0 1],[0 0],'k','LineWidth',2)
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.normalizedSpecificEnergy,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_bank_angle);
        end
        
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluatedBankAngle_v_normalizedSpecificEnergy_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end






%% Bank Reversal Trigger vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 122000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Reversal Trigger vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 2])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Bank Reversal Trigger (-)') % y-axis label
        set(gca,'YTick', 0:1:1);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger);
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankReversalTrigger_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Flight-Path Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        % for k = 30
        fig_num = p*100 + k*1 + 123000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Flight-Path Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-15 10])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Flight-Path Angle (deg)') % y-axis label
        set(gca,'YTick', -15:3:10);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        
        plot([0 60],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.flight_path_angle);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/flightPathAngle_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Flight-Path Angle Rate vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %  for k = 30
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 124000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Flight-Path Angle Rate vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(3*[-1 1])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Flight-Path Angle Rate (deg/s)') % y-axis label
        set(gca,'YTick', -3:1:3);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        plot([0 60],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.flight_path_angle_rate);
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/flightPathAngleRate_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end




%% Dynamic Pressure vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 30
        fig_num = p*100 + k*1 + 125000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Dynamic Pressure vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 10000]/1e3)
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Dynamic Pressure (kPa)') % y-axis label
        % set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        hold on
        grid on
        plot([0 60],[100 100]/1e3,'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.dynamicPressure/1e3);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/dynamicPressure_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Bending Moment vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 20
        fig_num = p*100 + k*1 + 126000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bending Moment vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 10])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Bending Moment (kPa-rad)') % y-axis label
        % set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        grid on
        plot([0 60],[5014 5014]/1e3,'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.bending_moment/1e3);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bendingMoment_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end




%% Pitch Moment Coefficient vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 127000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Cm vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-.001 .001])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('C_m (-)') % y-axis label
        set(gca,'YTick', -.001:0.0001:.001);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        plot([0 60],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.aero_moment_coefficient_C_m );
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Cm_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
    
end


%% BodyFlap Deflection Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        % for k = 1
        fig_num = p*100 + k*1 + 128000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('BodyFlap Deflection Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('BodyFlap Deflection Angle (deg)') % y-axis label
        set(gca,'YTick', -30:5:30);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        plot([0 60],[0 0],'k','LineWidth',2)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.bodyflap_deflection);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bodyflap_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Elevon Deflection Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 129000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Elevon Deflection Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-45 45])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Elevon Deflection Angle (deg)') % y-axis label
        set(gca,'YTick', -45:10:45);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        plot([0 60],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.elevon_deflection);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/elevon_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end




%% Total Body G-load vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 130000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Body G-load vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 2])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Total Body G-load (g)') % y-axis label
        set(gca,'YTick', 0:.2:2);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag );
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_body_g_load_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Height vs. Tauber Heat Flux at Leading Edge - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 131000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Tauber Heat Flux at Leading Edge - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 150])
        % xlim([0 8000])
        xlabel('Tauber Heat Flux at Leading Edge (W/m^2)') % x-axis label
        ylabel('Height (km)') % y-axis label
        % set(gca,'YTick', 0:10:150);
        % set(gca,'XTick', 0:500:8000);
        hold on
        
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_tauber_leadingedge,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_heatFluxTauber_leadingedge_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end

%% Tauber Heat Flux at Leading Edge vs. Airspeed - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 132000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Tauber Heat Flux at Leading Edge vs. Airspeed - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1000])
        xlim([0 8000])
        xlabel('Airspeed (m/s)') % x-axis label
        ylabel('Tauber Heat Flux at Leading Edge (W/m^2)') % y-axis label
        set(gca,'YTick', 0:100:1000);
        set(gca,'XTick', 0:500:8000);
        hold on
        
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.airspeed,...
                compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_tauber_leadingedge/1e3);
        end
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heatFluxTauber_leadingedge_v_V_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end

%% Tauber Heat Flux at Leading Edge vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 133000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Tauber Heat Rate at Leading Edge vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1000])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Tauber Heat Rate at Leading Edge (kW/m^2)') % y-axis label
        set(gca,'YTick', 0:100:1000);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_tauber_leadingedge/1e3 );
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heatRateTauber_leadingedge_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Chapman Heat Flux at Nose vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        %for k = 20
        fig_num = p*100 + k*1 + 134000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Chapman Heat Flux at Nose vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 500])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Chapman Heat Flux at Nose (kW/m^2)') % y-axis label
        set(gca,'YTick', 0:50:500);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_chapman_nose/1e3 );
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heatFluxChapman_nose_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Chapman Eq. Wall Temp. at Nose vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        %for k = 20
        fig_num = p*100 + k*1 + 135000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Chapman Eq. Wall Temp. at Nose vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 2000])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Chapman Eq. Wall Temp. at Nose (K)') % y-axis label
        set(gca,'YTick', 0:200:2000);
        set(gca,'XTick', 0:10:60);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.wall_temperature_chapman );
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/wallTempChapman_nose_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end




%% TUDAT Heat Rate at Nose vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 136000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('TUDAT Heat Rate at Nose vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1000])
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('TUDAT Heat Rate at Nose (kW/m^2)') % y-axis label
        set(gca,'YTick', 0:100:1000);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.heat_rate_TUDAT_nose/1e3 );
            
            
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                (compilation(p).evolutions(k).trajectories(ii).individual.heat_rate_TUDAT_nose - 0.5*(compilation(p).evolutions(k).trajectories(ii).individual.localDensity).*(compilation(p).evolutions(k).trajectories(ii).individual.airspeed).^3)/1e3 );
            
            
            
        end
        
        
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heatRateTUDAT_nose_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Height vs. Groundtrack - per Evolution
lat_f_deg = 38.9444444444444;
lon_f_deg = -77.4558333333;

lon_i_deg = -106.7;
lon_i_rad = deg2rad(-22.37);
lon_i_deg = -22.37;
lat_f_deg = 5;
lon_f_deg = -53;
validation = 1;


th = 0:pi/50:2*pi;
xunit = .75 * cos(th) + lon_f_deg;
yunit = .75 * sin(th) + lat_f_deg;
for p = 1:numel(compilation)
    
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 137000;
        figure(fig_num)
        hold on
        title(strcat('Trajectories - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('\tau (deg)') % x-axis label
        ylabel('\delta (deg)') % y-axis label
        zlabel('Height (km)') % z-axis label
        img = imread('img.jpg');
        imagesc([-180 180], [-90 90], (flipud(img)));
        set (gca,'Fontsize',20)
        set(gca,'XTick', -90:15:30);
        set(gca,'YTick', 0:15:90);
        set(gca,'ZTick', 0:15:150);
        xlim([-90 30])
        ylim([0 90])
        zlim([0 150])
        xlim([-82 -42])
        ylim([-8 12])
        
        
        for ii =1:numel(compilation(p).evolutions(k).trajectories)
            
            color_line3(...
                compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
                compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
                compilation(p).evolutions(k).trajectories(ii).individual.time_vector);
            stairs(...
                compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
                compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,'g')
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
            
        end
        
        plot(lon_f_deg,lat_f_deg,'MarkerSize',20)
        plot(lon_f_deg*[1 1],90*[-1 1],'k','LineWidth',2)
        plot(180*[-1 1],lat_f_deg*[1 1],'k','LineWidth',2)
        plot(xunit, yunit,'k','LineWidth',2);
        scatter(lon_f_deg,lat_f_deg,100,'r','x')
        axP = get(gca,'Position');
        %  legend(legendtext,'Location','southeastoutside')
        set(gca, 'Position', axP)
        view([13 49])
        grid on
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Height_v_Groundtrack_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set)',...
            '.png'),...
            'png');
    end
end




%% Height vs. Distance To Go vs. Mass - per Evolution
for p = 1:numel(compilation)
    
    
    % for k = 1:numel(compilation(p).evolutions)
    fig_num = p*100 + k*1 + 138000;
    figure(fig_num)
    hold on
    title(strcat('Height vs. Distance To Go vs. Mass - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
    set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
    xlabel('Distance To Go (deg)') % x-axis label
    ylabel('Mass (kg))') % y-axis label
    zlabel('Height (km)') % z-axis label
    % img = imread('img.jpg');
    %imagesc([-180 180], [-90 90], (flipud(img)));
    set (gca,'Fontsize',20)
    %set(gca,'XTick', -90:15:30);
    %set(gca,'YTick', 0:15:90);
    %set(gca,'ZTick', 0:15:150);
    %xlim([-90 30])
    %ylim([0 90])
    zlim([0 175])
    
    
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1        %     for k = 1:numel(compilation(p).evolutions)
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if sum(abs(compilation(p).evolutions(k).trajectories(ii).individual.heading_error) > 30) == 0
                if compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 30
                    if sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) > 2
                        
                        %color_line3(...
                        %    compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
                        %    compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,...
                        %    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
                        %    compilation(p).evolutions(k).trajectories(ii).individual.mass);
                        plot3(...
                            compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                            compilation(p).evolutions(k).trajectories(ii).individual.mass,...
                            compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
                        %                 scatter3(...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
                        %                 scatter3(...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
                    end
                end
            end
        end
    end
    
    axP = get(gca,'Position');
    %  legend(legendtext,'Location','southeastoutside')
    set(gca, 'Position', axP)
    view([13 49])
    grid on
    hold off
    saveas(...
        figure(fig_num),...
        strcat(...
        compilation(p).mainpath,...
        '/figures/Height_v_distanceToGo_v_mass_Set',...
        convertCharsToStrings(compilation(p).set)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end



%% Height vs. Distance To Go vs. Bank Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 139000;
        figure(fig_num)
        hold on
        title(strcat('Height vs. Distance To Go vs. Bank Angle - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Bank Angle (deg))') % y-axis label
        zlabel('Height (km)') % z-axis label
        % img = imread('img.jpg');
        %imagesc([-180 180], [-90 90], (flipud(img)));
        set (gca,'Fontsize',20)
        set(gca,'XTick', 0:10:60);
        % set(gca,'YTick', 0:15:90);
        set(gca,'ZTick', 0:25:175);
        xlim([0 60])
        %ylim([0 90])
        zlim([0 175])
        hold on
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if sum(abs(compilation(p).evolutions(k).trajectories(ii).individual.heading_error) > 30) == 0
                if compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 30
                    if sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) > 2
                        
                        plot3(...
                            compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                            compilation(p).evolutions(k).trajectories(ii).individual.bank_angle,...
                            compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
                        %plot3(...
                         %   compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                          %  compilation(p).evolutions(k).trajectories(ii).individual.skip_suppression_limit,...
                           % compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
                        
                        % stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                        %     compilation(p).evolutions(k).trajectories(ii).individual.bank_angle,'g');
                        % stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                        %     compilation(p).evolutions(k).trajectories(ii).individual.skip_suppression_limit,'k');
                        
                        % plot3(...
                        %     compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                        %     ones(numel(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go))*100,...
                        %     compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,'g')
                        
                        %                 scatter3(...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
                        %                 scatter3(...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
                        % end
                    end
                end
            end
        end
        
        axP = get(gca,'Position');
        %  legend(legendtext,'Location','southeastoutside')
        set(gca, 'Position', axP)
        view(3)
        %view([24 26])
        view([50 49])
        grid on
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Height_v_distanceToGo_v_bankAngle_Set',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set)',...
            '.png'),...
            'png');
        % close(fig_num);
        %  end
    end
    
end

%% Height vs. Distance To Go vs. Flight-Path Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 140000;
        figure(fig_num)
        hold on
        title(strcat('Height vs. Distance To Go vs. Flight-Path Angle - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Flight-Path Angle (deg))') % y-axis label
        zlabel('Height (km)') % z-axis label
        % img = imread('img.jpg');
        %imagesc([-180 180], [-90 90], (flipud(img)));
        set (gca,'Fontsize',20)
        %set(gca,'XTick', -90:15:30);
        set(gca,'YTick', -90:15:90);
        %set(gca,'ZTick', 0:15:150);
        xlim([0 60])
        %ylim([0 90])
        zlim([0 100])
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            
            %color_line3(...
            %    compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.mass);
            plot3(...
                compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.flight_path_angle,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
            
            
            %stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %   compilation(p).evolutions(k).trajectories(ii).individual.flight_path_angle,'g');
            
            
            % plot3(...
            %     compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %     ones(numel(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go))*100,...
            %     compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,'g')
            
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
            %end
            
            
            
        end
        
        axP = get(gca,'Position');
        %  legend(legendtext,'Location','southeastoutside')
        set(gca, 'Position', axP)
        view([24 26])
        grid on
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Height_v_distanceToGo_v_flightpathAngle_Set',...
            convertCharsToStrings(compilation(p).set)',...
            '.png'),...
            'png');
        % close(fig_num);
    end
end

%% Height vs. Distance To Go vs. Total Body G-load - per Evolution
for p = 1:numel(compilation)
    
    
    %for k = 1:numel(compilation(p).evolutions)
    for k = numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 141000;
        figure(fig_num)
        hold on
        title(strcat('Height vs. Distance To Go vs. Total Body G-load - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Total Body G-load (g)') % y-axis label
        zlabel('Height (km)') % z-axis label
        % img = imread('img.jpg');
        %imagesc([-180 180], [-90 90], (flipud(img)));
        set (gca,'Fontsize',20)
        set(gca,'XTick', 0:5:60);
        %set(gca,'YTick', 0:15:90);
        set(gca,'ZTick', 0:15:165);
        xlim([0 60])
        %ylim([0 90])
        zlim([0 165])
        
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            
            %color_line3(...
            %    compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.mass);
            plot3(...
                compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag);
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
        end
    end
    
    axP = get(gca,'Position');
    %  legend(legendtext,'Location','southeastoutside')
    set(gca, 'Position', axP)
    view([13 49])
    grid on
    hold off
    saveas(...
        figure(fig_num),...
        strcat(...
        compilation(p).mainpath,...
        '/figures/Height_v_distanceToGo_v_total_body_g_load_Set',...
        convertCharsToStrings(compilation(p).set)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end

%% Height vs. Distance To Go vs. Equilibrium Wall Temperature - per Evolution
for p = 1:numel(compilation)
    
    
    %for k = 1:numel(compilation(p).evolutions)
    for k = numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 142000;
        figure(fig_num)
        hold on
        title(strcat('Height vs. Distance To Go vs. Eq. Wall Temp. - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Eq. Wall Temp. (K)') % y-axis label
        zlabel('Height (km)') % z-axis label
        % img = imread('img.jpg');
        %imagesc([-180 180], [-90 90], (flipud(img)));
        set (gca,'Fontsize',20)
        set(gca,'XTick', 0:10:60);
        set(gca,'YTick', 0:100:1800);
        %set(gca,'ZTick', 0:15:150);
        %xlim([0 60])
        %ylim([0 90])
        zlim([0 175])
        
        grid on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            %      for ii = 1
            %color_line3(...
            %    compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.mass);
            plot3(...
                compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.wall_temperature_chapman,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.wall_temperature_chapman);
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
        end
    end
    
    axP = get(gca,'Position');
    %  legend(legendtext,'Location','southeastoutside')
    set(gca, 'Position', axP)
    view([13 49])
    grid on
    hold off
    saveas(...
        figure(fig_num),...
        strcat(...
        compilation(p).mainpath,...
        '/figures/Height_v_distanceToGo_v_total_body_g_load_Set',...
        convertCharsToStrings(compilation(p).set)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end



end
