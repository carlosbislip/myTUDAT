function [  ] = plotTimeHistories_AeroCoefficients( compilation )

%% Time History: Pitch Moment Coefficient - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 725000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Cm through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-.01 .01])
        max_tof = max([compilation(p).evolutions.max_tof]);
        %max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Time (s)') % x-axis label
        ylabel('C_m (-)') % y-axis label
        set(gca,'YTick', -.01:0.001:.01);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        plot([0 8000],(0)*[1 1],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.aero_moment_coefficient_C_m,'k','LineWidth',2);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.aero_moment_coefficient_C_m);
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Cm_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %        close(fig_num);
    end
end

%% Time History: Cm Increment due to Control Surfaces - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 724000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Cm Increment due to Control Surfaces through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-.1 .03])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Time (s)') % x-axis label
        ylabel('\Delta C_m (-)') % y-axis label
        set(gca,'YTick', -.1:0.01:.03);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        plot([0 8000],(0)*[1 1],'k','LineWidth',2)
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.increment_Cm_bodyflap_dif,'k','LineWidth',2);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.increment_Cm_bodyflap_dif);
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Cm_increment_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %        close(fig_num);
    end
end







end