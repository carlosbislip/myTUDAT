function [  ] = plotTrajectories( compilation, mainpath )


lat_f_deg = 38.9444444444444;
lon_f_deg = -77.4558333333;



% %% Termination State on 2D Map - per Set
% %close all
% for p = 1:numel(compilation)
%     fig_num = p + 100;
%     figure(fig_num)
%     hold on
%     title(strcat('Termination State - per Set',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
%     set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
%     xlabel('\tau (deg)') % x-axis label
%     ylabel('\delta (deg)') % y-axis label
%     img = imread('img.jpg');
%     imagesc([-180 180], [-90 90], (flipud(img)));
%     set (gca,'Fontsize',20)
%     set(gca,'YTick', -90:15:90);
%     set(gca,'XTick', -180:30:180);
%     xlim([-180 180])
%     ylim([-90 90])
%     for k = 1:numel(compilation(p).evolutions)
%         %         for ii = 1:numel(compilation(p).evolutions(k).trajectories)
%         %         terminal_lat(ii) =  compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle(end);
%         %                    lons(ii) =  compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle(end);
%         %
%         %
%         %         end
%         scatter(...
%             compilation(p).evolutions(k).individuals.lon_f_deg,...
%             compilation(p).evolutions(k).individuals.lat_f_deg,'x')
%         legendtext(k) = cellstr(strcat(strrep(convertCharsToStrings(compilation(p).set),'_',' '),' - Evolution:_{ } ', num2str(k-1)));
%     end
%
%     plot(lon_f_deg,lat_f_deg,'MarkerSize',20)
%     plot(lon_f_deg*[1 1],90*[-1 1],'k','LineWidth',2)
%     plot(180*[-1 1],lat_f_deg*[1 1],'k','LineWidth',2)

th = 0:pi/50:2*pi;
xunit = .75 * cos(th) + lon_f_deg;
yunit = .75 * sin(th) + lat_f_deg;
%     plot(xunit, yunit,'k','LineWidth',2);
%     scatter(lon_f_deg,lat_f_deg,100,'r','x')
%     axP = get(gca,'Position');
%     legend(legendtext,'Location','southeastoutside')
%     set(gca, 'Position', axP)
%
%     clear legendtext
%     %xlim([(lon_IAD_deg - 20) (lon_IAD_deg + 20)])
%     %ylim([(lat_IAD_deg - 10) (lat_IAD_deg + 10)])
%     hold off
%     saveas(...
%         figure(fig_num),...
%         strcat(...
%         mainpath,...
%         '/figures/Termination_State_2D_Map_per_Set',...
%         convertCharsToStrings(compilation(p).set),...
%         '.png'),...
%         'png');
%     close(fig_num);
%
% end
%
% %% Termination State on 2D Map - Zoom - per Set
%
% for p = 1:numel(compilation)
%     fig_num = p + 100;
%     figure(fig_num)
%     hold on
%     title(strcat('Termination State - Zoom - per Set ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
%     set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
%     xlabel('\tau (deg)') % x-axis label
%     ylabel('\delta (deg)') % y-axis label
%     img2 = imread('img.jpg');
%     imagesc([-180 180], [-90 90], (flipud(img2)));
%     set (gca,'Fontsize',20)
%     set(gca,'YTick', -90:1:90);
%     set(gca,'XTick', -180:1:180);
%     xlim([(lon_f_deg - 10) (lon_f_deg + 10)])
%     ylim([(lat_f_deg - 5) (lat_f_deg + 5)])
%
%     for k = 1:numel(compilation(p).evolutions)
%         scatter(...
%             compilation(p).evolutions(k).individuals.lon_f_deg,...
%             compilation(p).evolutions(k).individuals.lat_f_deg,'x')
%         legendtext(k) = cellstr(strcat(strrep(convertCharsToStrings(compilation(p).set),'_',' '),' - Evolution:_{ } ', num2str(k-1)));
%     end
%
%     plot(lon_f_deg*[1 1],90*[-1 1],'k','LineWidth',2)
%     plot(180*[-1 1],lat_f_deg*[1 1],'k','LineWidth',2)
%     plot(xunit, yunit,'k','LineWidth',2);
%     scatter(lon_f_deg,lat_f_deg,100,'r','x')
%     axP = get(gca,'Position');
%     legend(legendtext,'Location','southeastoutside')
%     set(gca, 'Position', axP)
%     clear legendtext
%
%     hold off
%
%     print(figure(fig_num),'MySavedPlot','-dpng')
%
%     saveas(...
%         figure(fig_num),...
%         strcat(...
%         mainpath,...
%         '/figures/Termination_State_2D_Map_Zoom_per_Set',...
%         convertCharsToStrings(compilation(p).set),...
%         '.png'),...
%         'png');
%     close(fig_num);
%
% end

%% Trajectories on 2D Map - per Set
%close all
for p = 1:numel(compilation)
    fig_num = p + 300;
    figure(fig_num)
    hold on
    title(strcat('Trajectories - per Set ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
    set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
    xlabel('\tau (deg)') % x-axis label
    ylabel('\delta (deg)') % y-axis label
    img = imread('img.jpg');
    imagesc([-180 180], [-90 90], (flipud(img)));
    set (gca,'Fontsize',20)
    set(gca,'YTick', -90:15:90);
    set(gca,'XTick', -180:30:180);
    xlim([-180 180])
    ylim([-90 90])
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(...
                compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
                compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle)
            legendtext(k) = cellstr(strcat(strrep(convertCharsToStrings(compilation(p).set),'_',' '),' - Evolution:_{ } ', num2str(k-1)));
        end
    end
    
    plot(lon_f_deg,lat_f_deg,'MarkerSize',20)
    plot(lon_f_deg*[1 1],90*[-1 1],'k','LineWidth',2)
    plot(180*[-1 1],lat_f_deg*[1 1],'k','LineWidth',2)
    plot(xunit, yunit,'k','LineWidth',2);
    scatter(lon_f_deg,lat_f_deg,100,'r','x')
    axP = get(gca,'Position');
    %  legend(legendtext,'Location','southeastoutside')
    set(gca, 'Position', axP)
    
    clear legendtext
    %xlim([(lon_IAD_deg - 20) (lon_IAD_deg + 20)])
    %ylim([(lat_IAD_deg - 10) (lat_IAD_deg + 10)])
    hold off
    saveas(...
        figure(fig_num),...
        strcat(...
        mainpath,...
        '/figures/Trajectories_2D_Map_Set',...
        convertCharsToStrings(compilation(p).set)',...
        '.png'),...
        'png');
   % close(fig_num);
    
end
%
% %% Trajectories on 2D Map - Zoom - per Set
% %close all
% for p = 1:numel(compilation)
%     fig_num = p + 310;
%     figure(fig_num)
%     hold on
%     title(strcat('Trajectories - Zoom - per Set ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
%     set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
%     xlabel('\tau (deg)') % x-axis label
%     ylabel('\delta (deg)') % y-axis label
%     img = imread('img.jpg');
%     imagesc([-180 180], [-90 90], (flipud(img)));
%     set (gca,'Fontsize',20)
%     set(gca,'YTick', -90:1:90);
%     set(gca,'XTick', -180:1:180);
%     xlim([(lon_f_deg - 10) (lon_f_deg + 10)])
%     ylim([(lat_f_deg - 5) (lat_f_deg + 5)])
%
%     if compilation(1).validation == 1
%         xlim([-72 -52])
%         ylim([-3 7])
%
%     end
%
%     for k = 1:numel(compilation(p).evolutions)
%         for ii = 1:numel(compilation(p).evolutions(k).trajectories)
%             plot(...
%                 compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
%                 compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle)
%             legendtext(k) = cellstr(strcat(strrep(convertCharsToStrings(compilation(p).set),'_',' '),' - Evolution:_{ } ', num2str(k-1)));
%         end
%     end
%
%     plot(lon_f_deg,lat_f_deg,'MarkerSize',20)
%     plot(lon_f_deg*[1 1],90*[-1 1],'k','LineWidth',2)
%     plot(180*[-1 1],lat_f_deg*[1 1],'k','LineWidth',2)
%     plot(xunit, yunit,'k','LineWidth',2);
%     scatter(lon_f_deg,lat_f_deg,100,'r','x')
%     axP = get(gca,'Position');
%     % legend(legendtext,'Location','southeastoutside')
%     set(gca, 'Position', axP)
%
%     clear legendtext
%     %xlim([(lon_IAD_deg - 20) (lon_IAD_deg + 20)])
%     %ylim([(lat_IAD_deg - 10) (lat_IAD_deg + 10)])
%     hold off
%     saveas...
%         (figure(fig_num),...
%         strcat(...
%         mainpath,...
%         '/figures/Trajectories_2D_Map_Zoom_Set',...
%         convertCharsToStrings(compilation(p).set),...
%         '.png'),...
%         'png');
%     close(fig_num);
%
% end
%
%
% %% Trajectory on 2D Map - Zoom - per Evolution
% %close all
% for p = 1:numel(compilation)
%
%     for k = 1:numel(compilation(p).evolutions)
%
%         fig_num = p*10 + 310 + k;
%         figure(fig_num)
%         hold on
%         title(strcat('Trajectories - Zoom - Evolution:_{ }',num2str(k - 1),' Set:_{ } ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
%         set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
%         xlabel('\tau (deg)') % x-axis label
%         ylabel('\delta (deg)') % y-axis label
%         img = imread('img.jpg');
%         imagesc([-180 180], [-90 90], (flipud(img)));
%         set (gca,'Fontsize',20)
%         set(gca,'YTick', -90:1:90);
%         set(gca,'XTick', -180:1:180);
%         xlim([(lon_f_deg - 10) (lon_f_deg + 10)])
%         ylim([(lat_f_deg - 5) (lat_f_deg + 5)])
%
%         if compilation(1).validation == 1
%
%             xlim([-72 -52])
%             ylim([-3 7])
%
%         end
%
%         for ii = 1:numel(compilation(p).evolutions(k).trajectories)
%
%             plot(...
%                 compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
%                 compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,'LineWidth',2)
%             %legendtext(k) = cellstr(strcat(strrep(convertCharsToStrings(compilation(p).set),'_',' '),' - Evolution:_{ } ', num2str(k-1)));
%
%         end
%
%         plot(lon_f_deg,lat_f_deg,'MarkerSize',20)
%         plot(lon_f_deg*[1 1],90*[-1 1],'k','LineWidth',2)
%         plot(180*[-1 1],lat_f_deg*[1 1],'k','LineWidth',2)
%         plot(xunit, yunit,'k','LineWidth',2);
%         scatter(lon_f_deg,lat_f_deg,100,'r','x')
%         axP = get(gca,'Position');
%         % legend(legendtext,'Location','southeastoutside')
%         set(gca, 'Position', axP)
%
%         clear legendtext
%         %xlim([(lon_IAD_deg - 20) (lon_IAD_deg + 20)])
%         %ylim([(lat_IAD_deg - 10) (lat_IAD_deg + 10)])
%         hold off
%         saveas(...
%             figure(fig_num),...
%             strcat(mainpath,...
%             '/figures/Trajectories_2D_Map_Zoom_Evolution_',...
%             num2str(k - 1),...
%             '_Set',...
%             convertCharsToStrings(compilation(p).set),...
%             '.png'),...
%             'png');
%         close(fig_num);
%
%     end
% end
%%
pp = 1;
'here';
%% Plot Trajectories on 3D Earth

for p = 1:numel(compilation)
    fig_num = p + 200 + pp;
    figure(fig_num)
    hold on
    title(strcat('Trajectory on spherical Earth - 1',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
    %set(figure(p + 200),'units','pixels','position',[0,0,1200,600])
    %   xlabel('\tau (deg)') % x-axis label
    %   ylabel('\delta (deg)') % y-axis label
    %   imagesc([-180 180], [-90 90], (flipud(img)));
    set (gca,'Fontsize',20)
    set(gca,'XTick',[], 'YTick', [])
    
    XTickLabel = [];
    YTickLabel = [];
    plotearth('Maptype','bluemarble',...
        'NEOMap',...
        strcat(mainpath,'additionalfunctions/PlotEarth/PlotEarth/BlueMarble.png'),...'/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/tudatBundle.git/tudatApplications/Space4ErrBody.git/Space4ErrBody_0/matlab/additionalfunctions/PlotEarth/PlotEarth/BlueMarble.png',...
        'SampleStep',1,'FullColor',true,'Shape','spheroid',...
        'viewpoint', [44 -5],fig_num,fig_num);
        axis('off');

    
    
    
        for i = 1:numel(compilation(p).evolutions)
            for k = 1:numel(compilation(p).evolutions(i).trajectories)
                color_line3(...
                    compilation(p).evolutions(i).trajectories(k).individual.x_R,...
                    compilation(p).evolutions(i).trajectories(k).individual.y_R,...
                    compilation(p).evolutions(i).trajectories(k).individual.z_R,...
                    compilation(p).evolutions(i).trajectories(k).individual.time_vector);
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
    view([44,-5])
    axis off
grid off
    ax = gca;               % get the current axis
    ax.Clipping = 'off';
    hold off
    saveas(figure(fig_num),strcat(mainpath,'/figures/Trajectories_3D_Map_1_Set',convertCharsToStrings(compilation(p).set),'.png'));
    
end
pp = pp + 1;
% %%
% for p = 1:numel(compilation)
% 
% figure(234234)
% hold on
% 
% 
%     
%         for i = 1:numel(compilation(p).evolutions)
%             for k = 1:numel(compilation(p).evolutions(i).trajectories)
%                 color_line3(...
%                     compilation(p).evolutions(i).trajectories(k).individual.x_R,...
%                     compilation(p).evolutions(i).trajectories(k).individual.y_R,...
%                     compilation(p).evolutions(i).trajectories(k).individual.z_R,...
%                     compilation(p).evolutions(i).trajectories(k).individual.time_vector);
%                 scatter3(...
%                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
%                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
%                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
%                 scatter3(...
%                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
%                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
%                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
%             end
%         end
%         
%         view(3)
% 
% end
% 
% 
% 
% 
% % %
% %
% %
% %
% %
% % %% Plot Trajectories on 3D Earth
% %
% % for p = 1:numel(compilation)
% %
% %     figure(p + 300)
% %     hold on
% %     title(strcat('Trajectory on spherical Earth - 2',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
% %     %set(figure(p + 200),'units','pixels','position',[0,0,1200,600])
% %     %   xlabel('\tau (deg)') % x-axis label
% %     %   ylabel('\delta (deg)') % y-axis label
% %     %   imagesc([-180 180], [-90 90], (flipud(img)));
% %     set (gca,'Fontsize',20)
% %     set(gca,'XTick',[], 'YTick', [])
% %
% %     axis('off');
% %     XTickLabel = [];
% %     YTickLabel = [];
% %     plotearth('Maptype','bluemarble','NEOMap',...
% %         '/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/code/PlotEarth/PlotEarth/BlueMarble.png',...
% %         'SampleStep',1,'FullColor',true,'Shape','spheroid',...
% %         'viewpoint', [-3 13],p + 300,p + 300);
% %
% %     for i = 1:numel(compilation(p).evolutions)
% %         for k = 1:numel(compilation(p).evolutions(i).trajectories)
% %             color_line3(...
% %                 compilation(p).evolutions(i).trajectories(k).individual.x_R,...
% %                 compilation(p).evolutions(i).trajectories(k).individual.y_R,...
% %                 compilation(p).evolutions(i).trajectories(k).individual.z_R,...
% %                 compilation(p).evolutions(i).trajectories(k).individual.time_vector);
% %             scatter3(...
% %                 compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
% %                 compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
% %                 compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
% %             scatter3(...
% %                 compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
% %                 compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
% %                 compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
% %         end
% %     end
% %     % view([12,-5])
% %     ax = gca;               % get the current axis
% %     ax.Clipping = 'off';
% %     hold off
% %     saveas(figure(fig_num),strcat(mainpath,'/figures/Trajectories_3D_Map_2_Set',convertCharsToStrings(compilation(p).set),'.png'));
% %
% % end
%
%
%
% %%
% % for p = 1:numel(compilation)
% %
% % %r = nan(1e6,26);
% % %t = nan(1e6,26);
% % figure(p + 400)
% % set(figure(p + 400),'units','pixels','position',[0,0,1200,600])
% % title('Trajectory in 3D - clipping disabled')
% % set (gca,'Fontsize',15)
% % xlabel('x (km) - Eath-Fixed') % x-axis label
% % ylabel('y (km) - Eath-Fixed') % y-axis label
% % zlabel('z (km) - Eath-Fixed') % z-axis label
% % xlim([4e3 5e3])
% % ylim([-6.5e3 0])
% % zlim([4.75e3 5.5e3])
% %
% % set(gca,'XTick', 4e3:.25e3:5e3);
% % set(gca,'YTick', -6.5e3:1e3:0e3);
% % set(gca,'zTick', 4.75e3:0.05e3:5.5e3);
% %
% % hold on
% % for k = 1:numel(compilation(p).evolutions(i).trajectories)
% %
% %
% %
% %
% %     plot3(compilation(p).evolutions(i).trajectories(k).individual.x_R/1e3,...
% %         compilation(p).evolutions(i).trajectories(k).individual.y_R/1e3,...
% %         compilation(p).evolutions(i).trajectories(k).individual.z_R/1e3);
% %
% %                 scatter3(compilation(p).evolutions(i).trajectories(k).individual.x_R(end)/1e3,...
% %                 compilation(p).evolutions(i).trajectories(k).individual.y_R(end)/1e3,...
% %                 compilation(p).evolutions(i).trajectories(k).individual.z_R(end)/1e3,'filled')
% %
% % end
% % view([44 -5])
% % ax = gca;               % get the current axis
% % ax.Clipping = 'off';
% % hold off
% %
% % saveas(gcf,strcat(mainpath,'/figures/Trajectory in 3D - 1 - ',strrep(convertCharsToStrings(compilation(p).set),'_',' '),'.png'));
% % end
%
% %% Grid Visualization: v_i, gamma_i, and chi_i
%
% % for p = 1:numel(compilation)
% %
% %     figure(p +500)
% %     hold on
% %     title(strcat('Inter-Continental Ballistic Tajectory - Cases tested -_{ } ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
% %     set(figure(p+500),'units','pixels','position',[0,0,1200,600])
% %     set (gca,'Fontsize',20)
% %        extents1 = 8e3;
% % extents2 = 90;
% % extents3 = 360;
% % %plot3(extents1*[0 1],[0 0],[0 0],'k')
% % %plot3([0 0],extents2*[0 1],[0 0],'k')
% % %plot3([0 0],[0 0],extents3*[0 1],'k')
% % xlim(extents1*[0 1])
% % ylim(extents2*[0 1])
% % %zlim(extents3*[0 1])
% % grid on
% % xlabel('v_i (deg)') % x-axis label
% % ylabel('\gamma_i (deg)') % y-axis label
% % zlabel('\chi_i (deg)') % y-axis label
% %      view([36,15])
% %
% %     for i = 1:numel(compilation(p).evolutions)
% %             scatter3(compilation(p).evolutions(i).v_i,...
% %                 compilation(p).evolutions(i).gamma_i,...
% %                 compilation(p).evolutions(i).chi_i);
% %
% %
% %     end
% %
% %     hold off
% %         saveas(gcf,strcat(mainpath,'/figures/Inter-Continental Ballistic Tajectory - Cases tested -_{ } ',strrep(convertCharsToStrings(compilation(p).set),'_',' '),'.png'));
% %
% % end
%
%






end