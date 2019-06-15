function [  ] = plotTimeHistories_AeroAngles( compilation )


%% Time History: Angle of Attack - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 1000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack  through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        max_tof = max([compilation(p).evolutions.max_tof]);
        %max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', 0:5:50);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
                
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack,'k','LineWidth',2);
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([10 45])
                set(gca,'YTick', 0:5:50);
                
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack);
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/AoA_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %   close(fig_num);
    end
end


%% Time History: Evaluated Angle of Attack - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 2000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Angle of Attack through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Evaluated Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', 0:5:50);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
                
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_angle_of_attack);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluated_AoA_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end

%% Time History: Commanded Angle of Attack - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Commanded Angle of Attack  through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', 0:5:50);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
           
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.commanded_angle_of_attack,'k','LineWidth',2);
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([10 45])
                set(gca,'YTick', 0:5:50);
                
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.commanded_angle_of_attack);
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/commanded_AoA_v_T_Evolution_',...
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
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 4000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-90 90])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.bank_angle,'k','LineWidth',2);
                plot([0 1400],[0 0],'k','LineWidth',2)
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([-80 100])
                set(gca,'YTick', -80:20:100);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bank_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end

%% Time History: Skip Suppresion Limit
for p = 1:numel(compilation)
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1
        fig_num = p*100 + 5000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Skip Suppression Limit through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(90*[-1 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Skip Suppression Limit (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],(0)*[1 1],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                abs(compilation(p).evolutions(k).trajectories(ii).individual.skip_suppression_limit));
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/skip_suppresion_limit_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %  close(fig_num);
    end
end





%% Time History: Evaluated Bank Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = 1
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 6000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Bank Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-90 90])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Evaluated Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_bank_angle);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluated_bank_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end

%% Time History: Commanded Bank Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = 1
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 7000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Commanded Bank Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-90 90])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.commanded_bank_angle);
        end
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/commanded_bank_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end

%% Time History: Bank Reversal Trigger - per Evolution
for p = 1:numel(compilation)
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 8000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Reversal Trigger - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1.5])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Bank Reversal Trigger (-)') % y-axis label
        set(gca,'YTick', 0:1:2);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        % for ii = numel(compilation(p).evolutions(k).trajectories):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger);
        end
        
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bank_reversal_trigger_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %  close(fig_num);
    end
end


%% Time History: Flight-Path Angle - per Evolution
for p = 1:numel(compilation)
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 9000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Flight-Path Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([-15 3])
        max_tof = max([compilation(p).evolutions.max_tof]);
        %max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Flight-Path Angle (deg)') % y-axis label
        %set(gca,'YTick', -15:3:3);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.flight_path_angle,'k','LineWidth',2);
                plot([0 1400],[0 0],'k','LineWidth',2)
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([-10 0])
                set(gca,'YTick', -10:1:0);
                
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.flight_path_angle);
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/flight_path_angle_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %  close(fig_num);
    end
end

%% Time History: Flight-Path Angle Rate - per Evolution
for p = 1:numel(compilation)
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 10000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Flight-Path Angle Rate through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-1 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        %max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Flight-Path Angle Rate (deg/s)') % y-axis label
        set(gca,'YTick', -1:.25:1);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.flight_path_angle_rate);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/flight_path_angle_rate_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %  close(fig_num);
    end
end





%% Time History: BodyFlap Deflection Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 11000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('BodyFlap Deflection Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-20 30])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof =1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('BodyFlap Deflection Angle (deg)') % y-axis label
        set(gca,'YTick', -20:5:30);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.bodyflap_deflection,'k','LineWidth',2);
                xlim([0 1400])
                        plot([0 1400],[0 0],'k','LineWidth',2)
                set(gca,'XTick', 0:200:1400);
                ylim([-20 20])
                set(gca,'YTick', -20:5:20);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.bodyflap_deflection);
            end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bodyflap_deflection_angle_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end



%% Time History: Elevon Deflection Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 12000 + k*1;
        
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Elevon Deflection Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-40 40])
        max_tof = max([compilation(p).evolutions.max_tof]);
        %max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Elevon Deflection Angle (deg)') % y-axis label
        set(gca,'YTick', -40:10:40);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.elevon_deflection,'k','LineWidth',2);
                xlim([0 1400])
                plot([0 1400],[0 0],'k','LineWidth',2)
                set(gca,'XTick', 0:200:1400);
                ylim([-40 40])
                set(gca,'YTick', -40:5:40);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.elevon_deflection);
            end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/elevon_deflection_angle_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end









end