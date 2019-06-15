function [  ] = plotIndividualsPerParameterAtNode( compilation )

%% Nodes - Ascent
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval - Ascent vs. Node - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 10])
        xlabel('Node (-)') % x-axis label
        ylabel('Node Interval - Ascent (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
   
        h = scatter(ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent1);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent2);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent3);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent4);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent5);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent6);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent7);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent8);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalAscent9);
        
        
        h = scatter(ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent1);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent2);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent3);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent4);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent5);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent6);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent7);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent8);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalAscent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/individualsPerNodeIntervalAtNodeAscent_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
         close(fig_num);
    end
end

%% Nodes - Descent
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Node Interval - Descent vs. Node - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        xlim([0 10])
        xlabel('Node (-)') % x-axis label
        ylabel('Node Interval - Descent (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
              h = scatter(ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent1);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent2);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent3);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent4);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent5);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent6);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent7);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent8);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.nodeIntervalDescent9);
        
        
        h = scatter(ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent1);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent2);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent3);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent4);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent5);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent6);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent7);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent8);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.nodeIntervalDescent9);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/individualsPerNodeIntervalAtNodeDescent_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
         close(fig_num);
    end
end

%% Angle of Attack - Ascent
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack - Ascent vs. Node - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([10 20])
        xlim([1 10])
        xlabel('Node (-)') % x-axis label
        ylabel('Angle of Attack - Ascent (deg)') % y-axis label
        set(gca,'YTick', 10:1:20);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        h = scatter(ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent1);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent2);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent3);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent4);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent5);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent6);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent7);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent8);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent9);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(10*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackAscent10);
        
        
        h = scatter(ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent1);
        
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent2);
        
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent3);
        
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent4);
        
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent5);
        
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent6);
        
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent7);
        
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent8);
        
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent9);
        
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(10*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackAscent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/individualsPerAngleOfAttackAtNodeAscent_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
         close(fig_num);
    end
end




%% Angle of Attack - Descent
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack - Descent vs. Node - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 50])
        xlim([1 10])
        xlabel('Node (-)') % x-axis label
        ylabel('Angle of Attack - Descent (deg)') % y-axis label
        set(gca,'YTick', 0:10:50);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
h = scatter(ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent1);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent2);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent3);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent4);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent5);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent6);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent7);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent8);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent9);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(10*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.angleOfAttackDescent10);
        
        
        h = scatter(ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent1);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent2);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent3);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent4);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent5);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent6);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent7);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent8);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent9);
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(10*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.angleOfAttackDescent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/individualsPerAngleOfAttackAtNodeDescent_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
         close(fig_num);
    end
end



%% Bank Angle - Descent
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 100000;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle - Descent vs. Node - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 85])
        xlim([1 10])
        xlabel('Node (-)') % x-axis label
        ylabel('Bank Angle - Descent (deg)') % y-axis label
        set(gca,'YTick', 0:10:85);
        %set(gca,'XTick', 0:2.5:30);
        hold on
        grid on
        
        h = scatter(ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent1);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent2);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent3);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent4);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent5);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent6);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent7);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent8);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent9);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(10*ones(compilation(p).evolutions(k).nonPrintedPopulationSize,1),...
            compilation(p).evolutions(k).nonPrintedPopulationDV.bankAngleDescent10);
        
        
        h = scatter(ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent1);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(2*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent2);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(3*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent3);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(4*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent4);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(5*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent5);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(6*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent6);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(7*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent7);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(8*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent8);
        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(9*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent9);
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter(10*ones(compilation(p).evolutions(k).printedPopulationSize,1),...
            compilation(p).evolutions(k).printedPopulationDV.bankAngleDescent10);
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/individualsPerBankAngleAtNodeDescent_v_Individual_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
         close(fig_num);
    end
end






end