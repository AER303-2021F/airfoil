function [h_fig] = Plotme(Plotmeinput)
%clear all
%clc
%format long


%% Define Variables........................................................

mode                = Plotmeinput.mode;
x                   = Plotmeinput.x;
y                   = Plotmeinput.y;
colorvalue          = Plotmeinput.colorvalue;
markerstyle         = Plotmeinput.markerstyle;
linestylevalue      = Plotmeinput.linestylevalue;
linewidthvalue      = Plotmeinput.linewidthvalue;
xlabelname          = Plotmeinput.xlabelname;
ylabelname          = Plotmeinput.ylabelname;
xtickvalue          = Plotmeinput.xtickvalue;
ytickvalue          = Plotmeinput.ytickvalue;
axisvalue           = Plotmeinput.axisvalue;
axisfontsize        = Plotmeinput.axisfontsize;
xtickdecimals       = Plotmeinput.xtickdecimals;
ytickdecimals       = Plotmeinput.ytickdecimals;
minorgridoption     = Plotmeinput.minorgridoption;
legendoption        = Plotmeinput.legendoption;
legendnames         = Plotmeinput.legend;
legendlocation      = Plotmeinput.legendlocation;
legendtextsize      = Plotmeinput.legendtextsize;
filename            = Plotmeinput.filename;
nVar                = Plotmeinput.nVar;
filelocation        = Plotmeinput.filelocation;
backlocation        = Plotmeinput.backlocation;
markersize          = Plotmeinput.markerSize;
markeredgecolor     = Plotmeinput.markeredgecolor;
markerfacecolor     = Plotmeinput.markerfacecolor;
textboxoption       = Plotmeinput.textboxoption;  
textboxtext         = Plotmeinput.textboxtext;     
textboxdim          = Plotmeinput.textboxdim;      
Positionvalue       = Plotmeinput.Positionvalue;
majorgridoption     = Plotmeinput.majorgridoption;

%% Color Pallate...........................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


colororder = [
    0.00  0.00  1.00
    0.00  0.50  0.00
    1.00  0.00  0.00
    0.00  0.75  0.75
    0.75  0.00  0.75
    0.75  0.75  0.00
    0.25  0.25  0.25
    0.75  0.25  0.25
    0.95  0.95  0.00
    0.25  0.25  0.75
    0.75  0.75  0.75
    0.00  1.00  0.00
    0.76  0.57  0.17
    0.54  0.63  0.22
    0.34  0.57  0.92
    1.00  0.10  0.60
    0.88  0.75  0.73
    0.10  0.49  0.47
    0.66  0.34  0.65
    0.99  0.41  0.23
    ];

%% Plot....................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if strcmp(mode,'Default') == 1;
    h_fig = figure;
    
    for i = 1 : nVar
        
        plot(x{1, i}, y{1, i}, 'Color', colorvalue{i}, ...
            'LineStyle', linestylevalue{i}, 'Linewidth', 2, ...
            'Marker', markerstyle{i}, 'MarkerFaceColor', colorvalue{i}, ...
            'MarkerSize', 3);
        hold on
        
    end
    
    hold off
    h_xlabel = xlabel(xlabelname);
    h_ylabel = ylabel(ylabelname);
    grid on;
    grid minor;
    
    set(gca,'FontSize',14)
    set(h_xlabel, 'Fontweight','bold', 'FontSize',30);
    set(h_ylabel, 'Fontweight','bold', 'FontSize',30);
    set(h_fig, 'Position', [100, 100, 1000, 700]);
    h_fig.PaperPositionMode = 'auto';
    
    ax = gca;
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.LineWidth = 1.5;
    ax.GridAlpha = 0.7;
    ax.MinorGridAlpha = 0.3;
    ax.XTick = xtickvalue;
    ax.YTick = ytickvalue;
    
    
    if strcmp(legendoption, 'On') == 1
        legendnames(legendnames, 'Location', 'best')
    end
    
    axis(axisvalue);
    print(filename,'-depsc','-r0')
    print(filename,'-dpng','-r0')
    clear h_xlabel h_ylabel ax
    
elseif strcmp(mode,'Advanced') == 1
    
    h_fig = figure;   
    for i = 1 : nVar
        
        plot(x{1, i}, y{1, i}, 'Color', colororder(colorvalue{i},:), ...
            'LineStyle', linestylevalue{i}, 'Linewidth', linewidthvalue, ...
            'Marker', markerstyle{i}, ...
            'MarkerFaceColor', colororder(markeredgecolor{i},:), ...
            'MarkerEdgeColor', colororder(markeredgecolor{i},:),...
            'MarkerSize', markersize{i});
        hold on
        
    end
    
    hold off
    h_xlabel = xlabel(xlabelname, 'Interpreter','latex');
    h_ylabel = ylabel(ylabelname, 'Interpreter','latex');
    if strcmp(majorgridoption, 'On') == 1
        grid on;
    end
   
    set(gca,'FontSize',legendtextsize)
    set(gca, 'FontName', 'Times New Roman')
    set(h_xlabel, 'Fontweight','bold', 'FontSize',axisfontsize);
    set(h_ylabel, 'Fontweight','bold', 'FontSize',axisfontsize);
    set(h_fig, 'Position', Positionvalue);
    h_fig.PaperPositionMode = 'auto';
    
    ax = gca;
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.LineWidth = 1.5;
    ax.GridAlpha = 0.7;
    
    if strcmp(minorgridoption, 'On') == 1
        grid minor;
        ax.MinorGridAlpha = 0.3;
    end
    
    ax.XTick = xtickvalue;
    ax.YTick = ytickvalue;
    
    if isstring(xtickdecimals) == 1
        if strcmp(xtickdecimals,'Default') == 1 
        else 
            disp('Plotme: Error input in xtickdecimals');   
        end
    else
        tixX=get(gca,'xtick')';
        set(gca,'xticklabel',num2str(tixX,['%.',num2str(xtickdecimals),'f']));
    end
    
    if isstring(ytickdecimals) == 1
        if strcmp(ytickdecimals,'Default') == 1
        else
            disp('Plotme: Error input in ytickdecimals');
        end
    else
        tixY=get(gca,'ytick')';
        set(gca,'yticklabel',num2str(tixY,['%.',num2str(ytickdecimals),'f']));
    end
    
    if strcmp(legendoption, 'On') == 1
        legend(legendnames, 'Location', legendlocation, 'Interpreter','latex')
    end
    
    axis(axisvalue);
    
    if strcmp(textboxoption, 'On') == 1
        annotation('textbox', textboxdim, 'String', textboxtext,...
            'FitBoxToText','on', 'Interpreter','latex',...
            'FontSize', legendtextsize, 'BackgroundColor', 'w');
    end
    
    
    if exist(filelocation, 'dir') == 7
        
        cd(filelocation)
        print(filename,'-depsc','-r0')
        print(filename,'-dpng','-r0')
        print(filename,'-dpdf','-r0')
        print(filename,'-dmeta','-r0')
        savefig([filename,'.fig']); 
        clear h_xlabel h_ylabel ax
        cd(backlocation)
        
    else
        
        mkdir(filelocation)
        cd(filelocation)
        print(filename,'-depsc','-r0')
        print(filename,'-dpng','-r0')
        print(filename,'-dpdf','-r0')
        print(filename,'-dmeta','-r0')
        savefig([filename,'.fig']); 
        clear h_xlabel h_ylabel ax
        cd(backlocation)
        
    end
    
else

    disp('Wrong Plotting Mode. Allowable options are Advanced and Default');
    
end



%% Demo....................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%{
clear all
clc
format long

XVAR  = -10 : 1 : 10;
YVAR = XVAR.^2;

Plotmeinput.x               = XVAR;
Plotmeinput.y               = YVAR;
Plotmeinput.colorvalue      = 'k';
Plotmeinput.xlabelname      = 'X';
Plotmeinput.ylabelname      = '\theta_0';
Plotmeinput.xtickvalue      = XVAR;
Plotmeinput.ytickvalue      = 0 : 10 : 100;
Plotmeinput.axisvalue       = [-10, 10, 0, 100];
Plotmeinput.foldername      = '';
Plotmeinput.filename        = 'testplot';
 
Plotme(Plotmeinput)
%}


%% Advanced Mode Plotting Lines............................................
%..........................................................................

%{
Plotmeinput.mode            = 'Advanced';
Plotmeinput.x               = {alfa_deg_interp{1,1}, alfa_deg_interp{1,2}};
Plotmeinput.y               = {imp_Cl{1,1}(:,4), imp_Cl{1,2}(:,4)};
Plotmeinput.colorvalue      = {1,7};
Plotmeinput.markerstyle     = {'none', 'none'};
Plotmeinput.markerSize      = {5,5};
Plotmeinput.markeredgecolor = {1,7};
Plotmeinput.markerfacecolor = {1,7};
Plotmeinput.linestylevalue  = {'-', '-.'};
Plotmeinput.linewidthvalue  = 3;
Plotmeinput.xlabelname      = '\alpha (deg)';
Plotmeinput.ylabelname      = 'C_l';
Plotmeinput.xtickvalue      = -10 : 5 : 20;
Plotmeinput.ytickvalue      = -0.5 : 0.5 : 1.5;
Plotmeinput.axisvalue       = [-10, 20, -0.5, 1.5];
Plotmeinput.Positionvalue   = [50, 50, 1000, 700];
Plotmeinput.axisfontsize    = 26;
Plotmeinput.xtickdecimals   = 0;
Plotmeinput.ytickdecimals   = 2;
Plotmeinput.majorgridoption = 'On';
Plotmeinput.minorgridoption = 'On';
Plotmeinput.legendoption    = 'On';
Plotmeinput.legend          = {'AS7004 Re = 100,000', 'AS7004 Re = 100,000 (Adjusted)'};
Plotmeinput.legendtextsize  = 24;
Plotmeinput.legendlocation  = 'best';
Plotmeinput.textboxoption   = 'Off';
Plotmeinput.textboxtext     = {'',''};
Plotmeinput.textboxdim      = []';                                          % Normalized units. Lower left corner (0, 0) and top right corner (1,1)
Plotmeinput.filename        = 'Comp_CLvsAlfa';
Plotmeinput.filelocation    = ' ';
Plotmeinput.backlocation    = ' ';
Plotmeinput.nVar            = 2;
Plotme(Plotmeinput)
%}



%% Credits................................................................

% Created by Suraj Bansal
% PhD Candidate
% Institute for Aerospace Studies
% University of Toronto
% Created: September 2016
