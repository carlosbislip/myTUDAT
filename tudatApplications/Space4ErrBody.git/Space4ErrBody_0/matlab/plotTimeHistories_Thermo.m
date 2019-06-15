function [  ] = plotTimeHistories_Thermo( compilation )





%% Time History: Chapman Heat Flux at Nose - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 654430 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Chapman Heat Flux at Nose through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
       %ylim([0 600])
        max_tof = max([compilation(p).evolutions.max_tof]);
        max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Chapman Heat Flux at Nose (kW/m^2)') % y-axis label
        set(gca,'YTick', 0:50:800);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on

           for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_chapman_nose/1e3,'k','LineWidth',2);
                            xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([0 600])
                set(gca,'YTick', 0:100:600);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_chapman_nose/1e3);
            end
        end
        
        plot([0 max_tof],530*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
compilation(p).mainpath,...
            '/figures/heatFluxChapman_nose_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
      %  close(fig_num);
    end
end

%% Time History: Chapman Wall Temp. at Nose - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 654430 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Chapman Wall Temp. at Nose through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 3000])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Chapman Wall Temp. at Nose (K)') % y-axis label
        set(gca,'YTick', 0:200:3000);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.wall_temperature_chapman,'k','LineWidth',2);
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([0 600])
                set(gca,'YTick', 0:100:600);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.wall_temperature_chapman);
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/wallTempChapman_nose_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %  close(fig_num);
    end
end



%% Time History: TUDAT Heat Rate at Nose - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 654430 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('TUDAT Heat Rate at Nose through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
       ylim([0 1000])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('TUDAT Heat Rate at Nose (kW/m)') % y-axis label
        set(gca,'YTick', 0:100:1000);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if  compilation(1).validation == 1
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.heat_rate_TUDAT_nose/1e3,'k','LineWidth',2);
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([0 600])
                set(gca,'YTick', 0:100:600);
            else
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.heat_rate_TUDAT_nose/1e3);
            end
        end

        %plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
compilation(p).mainpath,...
            '/figures/heatRateTUDAT_nose_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Time History: Tauber Heat Flux at Leading Edge - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 655430 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Tauber Heat Flux at Leading Edge through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1000])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Tauber Heat Flux at Leading Edge (kW/m)') % y-axis label
        set(gca,'YTick', 0:100:1000);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
for ii = 1:numel(compilation(p).evolutions(k).trajectories)   
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_tauber_leadingedge/1e3);
        end
        
        %plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
compilation(p).mainpath,...
            '/figures/heatFluxTauber_leadingedge_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end





end