function [  ] = plotParameterPerIndividual( compilation )

%{
populationSize = numel(compilation(1).evolutions(1).entirePopulation);
printedPopulationSize = numel(compilation(1).evolutions(1).printedPopulation);
nonPrintedPopulationSize = numel(compilation(1).evolutions(1).nonPrintedPopulation);

compilation(p).evolutions(k).populationSize = linspace(1,populationSize,populationSize);
filteredcompilation(p).evolutions(k).populationSize = compilation(p).evolutions(k).populationSize;%compilation(1).evolutions(1).filteredcompilation(p).evolutions(k).populationSize;


dif = compilation(p).evolutions(k).populationSize - filteredcompilation(p).evolutions(k).populationSize;



%AAAA = sum(~isnan(compilation(1).evolutions(1).filteredcompilation(p).evolutions(k).populationSize));
AAAA = sum(~isnan(filteredcompilation(p).evolutions(k).populationSize));
BBBB = sum(~isnan(dif));

%}

%% Node Interval Ascent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Ascent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Ascent 1 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalAscent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Ascent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Ascent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Ascent 2 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalAscent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Ascent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Ascent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Ascent 3 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalAscent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Ascent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Ascent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Ascent 4 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalAscent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Ascent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Ascent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Ascent 5 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalAscent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Ascent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Ascent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Ascent 6 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalAscent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Ascent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Ascent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Ascent 7 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalAscent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Ascent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Ascent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Ascent 8 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalAscent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Ascent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Ascent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Ascent 9 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalAscent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Angle of Attack Ascent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Ascent 1 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Angle of Attack Ascent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Ascent 2 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Angle of Attack Ascent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack 3 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Angle of Attack Ascent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Ascent 4 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Angle of Attack Ascent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Ascent 5 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Angle of Attack Ascent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Ascent 6 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Angle of Attack Ascent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Ascent 7 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Angle of Attack Ascent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Ascent 8 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Angle of Attack Ascent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Ascent 9 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Angle of Attack Ascent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Ascent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Ascent 10 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackAscent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 1 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 2 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 3 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 4 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 5 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 6 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 7 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 8 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 9 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Ascent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Ascent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Ascent 10 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleAscent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleAscent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleAscent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Thrust Elevation Angle Ascent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 1 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleAscent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Ascent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 2 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleAscent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Ascent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 3 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleAscent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Thrust Elevation Angle Ascent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 4 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleAscent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Ascent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 5 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleAscent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Thrust Elevation Angle Ascent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 6 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleAscent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Thrust Elevation Angle Ascent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 7 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngle7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Ascent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 8 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleAscent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Ascent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 9 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleAscent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Ascent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Ascent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Ascent 10 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleAscent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleAscent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleAscent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Ascent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 1 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Ascent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 2 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Ascent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 3 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Thrust Azimuth Angle Ascent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 4 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Ascent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 5 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Thrust Azimuth Angle Ascent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 6 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Thrust Azimuth Angle Ascent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 7 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Ascent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 8 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Ascent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 9 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Ascent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Ascent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Ascent 10 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleAscent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleAscent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleAscent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Throttle Setting Ascent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 1 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Throttle Setting Ascent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 2 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Throttle Setting Ascent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 3 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Throttle Setting Ascent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 4 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Throttle Setting Ascent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 5 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end





%% Throttle Setting Ascent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 6 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Throttle Setting Ascent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 7 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Throttle Setting Ascent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 8 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Throttle Setting Ascent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 9 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Throttle Setting Ascent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Ascent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Ascent 10 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingAscent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingAscent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingAscent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Initial Flight-Path Angle vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Initial Flight-Path Angle vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Initial Flight-Path Angle (deg)') % y-axis label
        set(gca,'YTick', -30:10:-30);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.initialFlightPathAngle);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.initialFlightPathAngle);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/initialFlightPathAngle_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Initial Launch Heading vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Initial Launch Heading vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 360])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Initial Launch Heading (deg)') % y-axis label
        set(gca,'YTick', 0:30:360);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.initialLaunchHeading);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.initialLaunchHeading);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/initialLaunchHeading_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Initial Velocity vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Initial Velocity vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Initial Velocity (m/s)') % y-axis label
        %set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.initialVelocity);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.initialVelocity);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/initialVelocity_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Maximum Velocity vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Maximum Velocity vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([3000 7000])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Maximum Velocity (m/s)') % y-axis label
        %set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.maximumVelocity);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.maximumVelocity);
        
        hold off
        %   set(gca,'xscale','log')
        %  set(gca,'yscale','log')
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/maximumVelocity_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Maximum Height vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Maximum Height vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Maximum Height (km)') % y-axis label
        %set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.maximumHeight/1e3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.maximumHeight/1e3);

        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/maximumHeight_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Additional Mass vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Additional Mass vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 100000]/1e3)
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Additional Mass (10^3kg)') % y-axis label
        %set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.additionalMass/1e3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.additionalMass/1e3);        
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/additionalMass_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Termination Distance Ratio vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Termination Distance Ratio vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Termination Distance Ratio (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.terminationDistanceRatio);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.terminationDistanceRatio);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/terminationDistanceRatio_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Node Interval Descent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Descent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Descent 1 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        %scatter(compilation(p).evolutions(k).populationSize,...
        %    compilation(p).evolutions(k).population.angleOfAttackDescent1);
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalDescent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Descent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Descent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Descent 2 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalDescent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Descent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Descent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Descent 3 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalDescent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Descent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Descent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Descent 4 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalDescent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Descent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Descent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Descent 5 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalDescent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Descent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Descent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Descent 6 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalDescent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Descent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Descent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Descent 7 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalDescent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Descent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Descent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Descent 8 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalDescent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Node Interval Descent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval Descent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Node Interval Descent 9 (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/nodeIntervalDescent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Angle of Attack Descent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Descent 1 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Angle of Attack Descent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Descent 2 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Angle of Attack Descent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack 3 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Angle of Attack Descent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Descent 4 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Angle of Attack Descent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Descent 5 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Angle of Attack Descent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Descent 6 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Angle of Attack Descent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Descent 7 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Angle of Attack Descent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Descent 8 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Angle of Attack Descent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Descent 9 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Angle of Attack Descent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack Descent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Angle of Attack Descent 10 (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleofAttackDescent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 1 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 2 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 3 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 4 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 5 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 6 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 7 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 8 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 9 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Bank Angle Descent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle Descent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 90])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Bank Angle Descent 10 (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngleDescent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Thrust Elevation Angle Descent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 1 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleDescent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Descent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 2 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleDescent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Descent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 3 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleDescent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Thrust Elevation Angle Descent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 4 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleDescent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Descent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 5 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleDescent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Thrust Elevation Angle Descent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 6 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleDescent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Thrust Elevation Angle Descent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 7 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngle7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Descent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 8 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleDescent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Descent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 9 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleDescent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Elevation Angle Descent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Elevation Angle Descent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Elevation Angle Descent 10 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustElevationAngleDescent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustElevationAngleDescent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustElevationAngleDescent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Descent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 1 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Descent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 2 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Descent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 3 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Thrust Azimuth Angle Descent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 4 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Descent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 5 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Thrust Azimuth Angle Descent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 6 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Thrust Azimuth Angle Descent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 7 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Descent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 8 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Descent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 9 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Thrust Azimuth Angle Descent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Azimuth Angle Descent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-25 25])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Thrust Azimuth Angle Descent 10 (deg)') % y-axis label
        set(gca,'YTick',-25:5:25);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        plot([0 compilation(p).evolutions(k).populationSize(end)],(0)*[1 1],'k','LineWidth',2)
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.thrustAzimuthAngleDescent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.thrustAzimuthAngleDescent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustAzimuthAngleDescent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Throttle Setting Descent 1 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 1 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 1 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent1);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent1);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent1_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Throttle Setting Descent 2 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 2 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 2 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent2);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent2);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent2_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Throttle Setting Descent 3 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 3 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 3 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent3);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent3);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent3_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Throttle Setting Descent 4 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 4 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 4 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent4);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent4);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent4_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Throttle Setting Descent 5 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 5 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 5 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent5);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent5);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent5_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end





%% Throttle Setting Descent 6 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 6 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 6 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent6);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent6);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent6_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Throttle Setting Descent 7 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 7 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 7 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent7);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent7);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent7_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Throttle Setting Descent 8 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 8 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 8 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent8);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent8);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent8_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end




%% Throttle Setting Descent 9 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 9 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 9 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent9);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent9_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Throttle Setting Descent 10 vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting Descent 10 vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Throttle Setting Descent 10 (-)') % y-axis label
        set(gca,'YTick',0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.throttleSettingDescent10);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.throttleSettingDescent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttleSettingDescent10_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end






%% Final Velocity vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Final Velocity vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Final Velocity (m/s)') % y-axis label
        %set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.finalVelocity);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.finalVelocity);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/finalVelocity_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Skip Suppression Trigger Time vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Skip Suppression Trigger Time vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Time (s)') % y-axis label
        %set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.finalVelocity);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.skipSuppressionTriggerTime);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/skipSuppressionTriggerTime_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Maximum Mechanical Load vs. Individual
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Maximum Mechanical Load vs. Individual - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 50])
        xlim([0 compilation(p).evolutions(k).populationSize])
        xlabel('Individual (-)') % x-axis label
        ylabel('Maximum Mechanical Load (g)') % y-axis label
        %set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        scatter(compilation(p).evolutions(k).nonPrintedIndices,...
            compilation(p).evolutions(k).nonPrintedPopulationDV.maximumMechanicalLoad);
        scatter(compilation(p).evolutions(k).printedIndices,...
            compilation(p).evolutions(k).printedPopulationDV.maximumMechanicalLoad  );
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/maximumMechanicalLoad_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



end


