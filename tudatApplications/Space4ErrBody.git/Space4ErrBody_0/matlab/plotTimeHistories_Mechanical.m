function [  ] = plotTimeHistories_Mechanical( compilation )


%% Time History: Total Body G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + 3957000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Body G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 3])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Total Body G-load (g)') % y-axis label
        %set(gca,'YTick', 0:.2:5);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag);
        end
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_body_g_load_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %        close(fig_num);
    end
end



%% Time History: Total Body Z-Component G-load - per Evolution
for p = 1:numel(compilation)
    
    %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1
        fig_num = p*200 + 3457000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Body Z-Component G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-5 10])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Total Body Z-Component G-load (g)') % y-axis label
        set(gca,'YTick', -5:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        % for ii = (numel(compilation(p).evolutions(k).trajectories)):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            %if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_z);
            % end
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_body_g_load_z_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end


%% Time History: Total Passenger Z-Component G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*200 + 3458000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Passenger Z-Component G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-5 10])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Total Passenger Z-Component G-load (g)') % y-axis label
        set(gca,'YTick', -5:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        grid on
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.passenger_fixed_total_g_load_z);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_passenger_g_load_z_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end



%% Time History: Aero G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3457000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Aero G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 2])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Aero G-load (g)') % y-axis label
        set(gca,'YTick', 0:.1:2);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.total_aero_g_load,'k','LineWidth',2);
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([0 2])
                set(gca,'YTick', 0:.2:2);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.total_aero_g_load);
            end
        end
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/aero_g_load_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end


%% Time History: Lift Force Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3485100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Lift Force Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 15])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Lift Force Magnitude (10^3 kN)') % y-axis label
        set(gca,'YTick', 0:1:15);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.currentLiftForce/1e6);
        end
        %legend('x-dir','y-dir','z-dir')
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/currentLiftForce_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %   close(fig_num);
    end
end


%% Time History: Drag Force Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3475100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Drag Force Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 15])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Drag Force Magnitude (10^3 kN)') % y-axis label
        set(gca,'YTick', 0:1:15);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.currentDragForce/1e6);
        end
        %legend('x-dir','y-dir','z-dir')
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/currentDragForce_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %   close(fig_num);
    end
end


%% Time History: Thrust Acc. Components - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3465000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Components through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 5])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Thrust Acc. Components (m/s^2)') % y-axis label
        %set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 30 )
                % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
                plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_x,'k');
                plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_y,'r');
                plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_z,'b');
            end
        end
        legend('x-dir','y-dir','z-dir')
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrust_acc_components_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Time History: Thrust Acc. Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3465100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Acc. Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 5])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Thrust Acc. Magnitude (m/s^2)') % y-axis label
        %set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 30 )
                plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_M);
            end
        end
        %legend('x-dir','y-dir','z-dir')
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Thrust_acc_magnitude_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %   close(fig_num);
    end
end









%% Time History: Thrust Force Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3465100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Force Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 15])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Thrust Force Magnitude (10^3 kN)') % y-axis label
        set(gca,'YTick', 0:1:15);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.thrustMagnitude/1e6);
        end
        %legend('x-dir','y-dir','z-dir')
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustMagnitude_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %   close(fig_num);
    end
end

%% Time History: Thrust g-Load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3465100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust g-Load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 10])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Thrust g-Load (g)') % y-axis label
        set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.totalThrustGLoad);
        end
        %legend('x-dir','y-dir','z-dir')
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustMagnitude_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %   close(fig_num);
    end
end



%% Time History: Gravity Force Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3466300 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Gravity Force Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 10])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Gravity Force Magnitude (N)') % y-axis label
        set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        % for ii = (numel(compilation(p).evolutions(k).trajectories)):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            %  if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                ((compilation(p).evolutions(k).trajectories(ii).individual.local_gravity_1).^2 + (compilation(p).evolutions(k).trajectories(ii).individual.local_gravity_2).^2).^(1/2));
            %  end
        end
        %legend('x-dir','y-dir','z-dir')
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/gravityMagnitude_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %   close(fig_num);
    end
end



%% Time History: Dynamic Pressure - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        % for k = 20
        fig_num = p*100 + 654180 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Dynamic Pressure through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 10])
        max_tof = max([compilation(p).evolutions.max_tof]);
        %max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Dynamic Pressure (kPa)') % y-axis label
        set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.dynamicPressure/1e3,'k','LineWidth',2);
                           xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                %ylim([0 1])
                %set(gca,'YTick', 0:.2:2);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.dynamicPressure/1e3);
            end
        end
        
        %plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/dynamic_pressure_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end

%% Time History: Bending Moment (Q-alpha) - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        % for k = 20
        fig_num = p*100 + 674180 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bending Moment (Q-alpha) through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 6])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Bending Moment, aka (Q-alpha) (kPa-rad)') % y-axis label
        %set(gca,'YTick', -90:30:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        plot([0 max_tof],(5014)*[1 1]/1e3,'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.bending_moment/1e3,'k','LineWidth',2);
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bending_moment_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end



end