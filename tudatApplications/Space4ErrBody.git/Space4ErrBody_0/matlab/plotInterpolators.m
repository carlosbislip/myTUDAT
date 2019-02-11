function [  ] = plotInterpolators( compilation, mainpath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

lat_f_deg = 38.9444444444444;
lon_f_deg = -77.4558333333;


% Coordinates for Validation case: Re-entry towards Kourou
if compilation(1).validation == 1
    
    lon_i_deg = -106.7;
    lon_i_rad = deg2rad(lon_i_deg);
    lat_f_deg = 5;
    lon_f_deg = -53;
    validation = 1;
    
end


%% Interpolators: Angle of Attack - Ascent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Interpolated AoA for Ascent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 30])
        max_interp_E_mapped_Ascent = max([compilation(p).evolutions.max_interp_E_mapped_Ascent]);
        xlim([0 max_interp_E_mapped_Ascent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Angle of attack (deg)') % y-axis label
        set(gca,'YTick', 0:1:30);
        set(gca,'XTick', 0:.2:max_interp_E_mapped_Ascent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_angle_of_attack_Ascent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_angle_of_attack_Ascent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_AoA_Ascent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Interpolators: Angle of Attack - Descent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Interpolated AoA for Descent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        max_interp_E_mapped_Descent = max([compilation(p).evolutions.max_interp_E_mapped_Descent]);
        xlim([0 max_interp_E_mapped_Descent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Angle of attack (deg)') % y-axis label
        set(gca,'YTick', 0:5:50);
          set(gca,'XTick', 0:max_interp_E_mapped_Descent/10:max_interp_E_mapped_Descent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_angle_of_attack_Descent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_angle_of_attack_Descent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_AoA_Descent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Interpolators: Bank Angle - Ascent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3439000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Interpolated Bank Angle for Ascent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 5])
        max_interp_E_mapped_Ascent = max([compilation(p).evolutions.max_interp_E_mapped_Ascent]);
        xlim([0 max_interp_E_mapped_Ascent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Bank angle (deg)') % y-axis label
        set(gca,'YTick', 0:1:5);
        set(gca,'XTick', 0:max_interp_E_mapped_Ascent/10:max_interp_E_mapped_Ascent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_bank_angle_Ascent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_bank_angle_Ascent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_bank_Ascent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Interpolators: Bank Angle - Descent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3439000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Interpolated Bank Angle for Descent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        max_interp_E_mapped_Descent = max([compilation(p).evolutions.max_interp_E_mapped_Descent]);
        %xlim([0 max_interp_E_mapped_Descent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Bank angle (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        set(gca,'XTick', 0:max_interp_E_mapped_Descent/10:max_interp_E_mapped_Descent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_bank_angle_Descent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_bank_angle_Descent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_bank_Descent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Interpolators: Thrust Elevation Angle - Ascent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Interpolated eps_T for Ascent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        max_interp_E_mapped_Ascent = max([compilation(p).evolutions.max_interp_E_mapped_Ascent]);
        xlim([0 max_interp_E_mapped_Ascent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Thrust Elevation Angle (deg)') % y-axis label
        % set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:.2:max_interp_E_mapped_Ascent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_thrust_elevation_angle_Ascent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_thrust_elevation_angle_Ascent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_eps_T_Ascent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Interpolators: Thrust Elevation Angle - Descent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Interpolated eps_T for Descent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        %max_interp_E_mapped_Descent = max(abs([compilation(p).evolutions.max_interp_E_mapped_Descent]));
        %xlim([0 max_interp_E_mapped_Descent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Thrust Elevation Angle (deg)') % y-axis label
        % set(gca,'YTick', 0:1:10);
       % set(gca,'XTick', 0:.2:max_interp_E_mapped_Descent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_thrust_elevation_angle_Descent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_thrust_elevation_angle_Descent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_eps_T_Descent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end





%% Interpolators: Thrust Azimuth Angle - Ascent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459300 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Interpolated phi_T for Ascent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        %max_interp_E_mapped_Descent = max(abs([compilation(p).evolutions.max_interp_E_mapped_Descent]));
        %xlim([0 max_interp_E_mapped_Descent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Thrust Azimuth Angle (deg)') % y-axis label
        % set(gca,'YTick', 0:1:10);
       % set(gca,'XTick', 0:.2:max_interp_E_mapped_Descent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_thrust_azimuth_angle_Ascent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_thrust_azimuth_angle_Ascent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_phi_T_Ascent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Interpolators: Thrust Azimuth Angle - Descent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459300 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Interpolated phi_T for Descent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        %max_interp_E_mapped_Descent = max(abs([compilation(p).evolutions.max_interp_E_mapped_Descent]));
        %xlim([0 max_interp_E_mapped_Descent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Thrust Azimuth Angle (deg)') % y-axis label
        % set(gca,'YTick', 0:1:10);
       % set(gca,'XTick', 0:.2:max_interp_E_mapped_Descent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_thrust_azimuth_angle_Descent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_thrust_azimuth_angle_Descent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_phi_T_Descent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Interpolators: Throttle Setting - Ascent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting for Ascent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-1 2])
        % max_tof = max([compilation(p).evolutions.max_tof]);
%        max_interp_E_mapped_Ascent = max([compilation(p).evolutions.max_interp_E_mapped]);
       % xlim([0 max_interp_E_mapped_Ascent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Throttle Setting (-)') % y-axis label
        % set(gca,'YTick', 0:1:10);
        %set(gca,'XTick', 0:.2:max_interp_E_mapped_Ascent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_throttle_setting_Ascent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_throttle_setting_Ascent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_throttle_Ascent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Interpolators: Throttle Setting - Descent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting for Descent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-1 2])
        % max_tof = max([compilation(p).evolutions.max_tof]);
%        max_interp_E_mapped_Ascent = max([compilation(p).evolutions.max_interp_E_mapped]);
       % xlim([0 max_interp_E_mapped_Ascent])
        xlabel('Mapped Energy: E_{mapped}') % x-axis label
        ylabel('Throttle Setting (-)') % y-axis label
        % set(gca,'YTick', 0:1:10);
        %set(gca,'XTick', 0:.2:max_interp_E_mapped_Ascent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.interp_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.interp_throttle_setting_Descent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_throttle_setting_Descent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/interp_throttle_Descent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Interpolators: Node Location - Ascent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Location for Ascent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        % max_tof = max([compilation(p).evolutions.max_tof]);
%        max_interp_E_mapped_Ascent = max([compilation(p).evolutions.max_interp_E_mapped]);
        xlim([0 1])
        ylabel('Mapped Energy: E_{mapped}') % x-axis label
        xlabel('Node Location (-)') % y-axis label
        % set(gca,'YTick', 0:1:10);
        %set(gca,'XTick', 0:.2:max_interp_E_mapped_Ascent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.DV_node_location_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Ascent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_node_location_Ascent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Ascent);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/node_location_Ascent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Interpolators: Node Location - Descent
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3459000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Location for Descent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
       % ylim([0 1])
        % max_tof = max([compilation(p).evolutions.max_tof]);
%        max_interp_E_mapped_Ascent = max([compilation(p).evolutions.max_interp_E_mapped]);
        %xlim([0 1])
        ylabel('Mapped Energy: E_{mapped}') % x-axis label
        xlabel('Node Location (-)') % y-axis label
        % set(gca,'YTick', 0:1:10);
        %set(gca,'XTick', 0:.2:max_interp_E_mapped_Ascent);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.DV_node_location_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Descent);
            scatter(compilation(p).evolutions(k).trajectories(ii).individual.DV_node_location_Descent,...
                compilation(p).evolutions(k).trajectories(ii).individual.DV_E_mapped_Descent);
        end
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/node_location_Descent_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


end
