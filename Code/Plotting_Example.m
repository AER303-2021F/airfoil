clear all
clc
format long
close all
warning('off','all')


% Let's generate 4 different functions of x................................

x{1} = linspace(0,1,100);
x{2} = linspace(0,pi, 10);

y{1} = x{2};
y{2} = x{1}.^2;
y{3} = sin(x{2});
y{4} = log(x{1});

% Advantage of using a cell is that the cell contents can hold arrays of
% different sizes.

% Let's plot the first function only.......................................

% Creating a structure "Plotmeinput" with different fields that gets
% passsed on to the function Plotme.m. These fields can be tweaked to suit
% your plotting preferences. This is a good method to make your
% plotting more efficient (and make them look professional). 

Plotmeinput.mode            = 'Advanced';                                   % Specifies the plotting mode. Leave it to "Advanced".
Plotmeinput.x               = {x{1}};                                       % If you want to plot n data sets, mention the x-axis variables for those data sets on this line
Plotmeinput.y               = {y{4}};                                       % If you want to plot n data sets, mention the x-axis variables for those data sets on this line
Plotmeinput.colorvalue      = {1};                                          % Specify color for each of the n data sets plot lines. For n data sets, you will need n color values. The n color values can all be same or different or some combination of that. Will need to experiment with numbers to see which colors you prefer using.
Plotmeinput.markerstyle     = {'None'};                                        % Specify the marker type here. See matlab documentation for supported markerstyles. For no markers, use 'None'.
Plotmeinput.markerSize      = {5};                                          % Mention the marker size here. In case of no markers, the values mentioned here are not used, but are still needed for the code to work
Plotmeinput.markeredgecolor = {7};                                          % Mention the marker edge color here. In case of no markers, the values mentioned here are not used, but are still needed for the code to work
Plotmeinput.markerfacecolor = {7};                                          % Mention the marker face color here. In case of no markers, the values mentioned here are not used, but are still needed for the code to work
Plotmeinput.linestylevalue  = {'-'};                                        % Specify the line type here. See matlab documentation for supported linestyles. In case you only want to point markers without lines, mention 'None'
Plotmeinput.linewidthvalue  = 3;                                            % Line thickness can be specified here. 
Plotmeinput.xlabelname      = '$x$';                                        % X-axis labels go here. Supports LaTeX inputs. Use $$ signs for latex math commands
Plotmeinput.ylabelname      = '$y(x)$';                                     % Y-axis labels go here. Supports LaTeX inputs. Use $$ signs for latex math commands                            
Plotmeinput.xtickvalue      = 0:0.2:1;                                    % Specify the x-tick values. 
Plotmeinput.ytickvalue      = -5 : 0.5: 0;                                 % Specify the y-tick values. 
Plotmeinput.axisvalue       = [0, 1, -5, 0];                             % Specify the axis limits
Plotmeinput.Positionvalue   = [50, 50, 1000, 700];                          % Changes the aspect ratio of plots. This might be advanced for some. If so, leave it at default. 
Plotmeinput.axisfontsize    = 26;                                           % Changes the font size of the numbers and labels on the axes. 
Plotmeinput.xtickdecimals   = 4;                                            % Specify how many decimal places you want on your x-axis numbers
Plotmeinput.ytickdecimals   = 2;                                            % Specify how many decimal places you want on your y-axis numbers
Plotmeinput.majorgridoption = 'On';                                         % Major grid on -> 'On' (case sensitive). Major grid off -> Any keyword except 'On'. Use 'Off' to avoid confusion later. 
Plotmeinput.minorgridoption = 'Off';                                        % Minor grid on -> 'On' (case sensitive). Minor grid off -> Any keyword except 'On'. Use 'Off' to avoid confusion later. 
Plotmeinput.legendoption    = 'Off';                                         % Want a legend? -> 'On' (case sensitive). Legend off -> Any keyword except 'On'. Use 'Off' to avoid confusion later. 
Plotmeinput.legend          = {'Linear', 'Quadratic',  'Sine',  'Logarithm'};% For n number of data sets, n number of legend entries are required. In this tutorial, n = 1.
Plotmeinput.legendtextsize  = 24;                                           % Specify the legend font size here.                                      
Plotmeinput.legendlocation  = 'best';                                       % Specify legend location. "Best" works fine most times. More more legend location options, consult MATLAB documentation. 
Plotmeinput.textboxoption   = 'Off';                                        % Want to put a text box? -> 'On'. You know by now what to do to leave it at off. Might be advanced. Don't use it if not familiar. 
Plotmeinput.textboxtext     = {'',''};                                      % Specify the text-box text here.
Plotmeinput.textboxdim      = []';                                          % Normalized units. Lower left corner (0, 0) and top right corner (1,1)
Plotmeinput.filename        = 'Example_Linear';                             % Specify the file name for the figure. 
Plotmeinput.filelocation    = ' ';                                          % To save the file in the same folder as this MATLAB script use ' ' (needs to have one space between two apostrophes). To save it in a different folder, specify the FilePath. If the folders don't exist, this script will automatically create them.
Plotmeinput.backlocation    = ' ';                                          % Need to specify the return path to the folder from which this script is running to avoid errors. 
Plotmeinput.nVar            = 1;                                            % Specify the number 'n' here. In this tutorial n = 1, where n is the number of datasets you want to co-plot on a single figure.
Plotme(Plotmeinput);                    


% Let's co-plot the four functions on the same plot for comparison.........

Plotmeinput.mode            = 'Advanced';                                   % Specifies the plotting mode. Leave it to "Advanced".
Plotmeinput.x               = {x{2}, x{1}, x{2}, x{1}};                     % If you want to plot n data sets, mention the x-axis variables for those data sets on this line
Plotmeinput.y               = {y{1}, y{2}, y{3}, y{4}};                     % If you want to plot n data sets, mention the x-axis variables for those data sets on this line
Plotmeinput.colorvalue      = {7, 1, 3, 5, 9};                              % Specify color for each of the n data sets plot lines. For n data sets, you will need n color values. The n color values can all be same or different or some combination of that. Will need to experiment with numbers to see which colors you prefer using.
Plotmeinput.markerstyle     = {'^', 'None', 'o', 'None'};                   % Specify the marker type here. See matlab documentation for supported markerstyles. For no markers, use 'None'.
Plotmeinput.markerSize      = {5, 5, 5, 5};                                 % Mention the marker size here. In case of no markers, the values mentioned here are not used, but are still needed for the code to work
Plotmeinput.markeredgecolor = {7, 1, 3, 5, 9};                              % Mention the marker edge color here. In case of no markers, the values mentioned here are not used, but are still needed for the code to work
Plotmeinput.markerfacecolor = {7, 1, 3, 5, 9};                              % Mention the marker face color here. In case of no markers, the values mentioned here are not used, but are still needed for the code to work
Plotmeinput.linestylevalue  = {'-', '-.', 'None', '--'};                    % Specify the line type here. See matlab documentation for supported linestyles. In case you only want to point markers without lines, mention 'None'
Plotmeinput.linewidthvalue  = 3;                                            % Line thickness can be specified here. 
Plotmeinput.xlabelname      = '$x$';                                        % X-axis labels go here. Supports LaTeX inputs. Use $$ signs for latex math commands
Plotmeinput.ylabelname      = '$y(x)$';                                     % Y-axis labels go here. Supports LaTeX inputs. Use $$ signs for latex math commands                            
Plotmeinput.xtickvalue      = 0:0.5:4;                                      % Specify the x-tick values. 
Plotmeinput.ytickvalue      = -6 : 2: 6;                                    % Specify the y-tick values. 
Plotmeinput.axisvalue       = [0, 4, -6, 6];                                % Specify the axis limits
Plotmeinput.Positionvalue   = [50, 50, 1000, 700];                          % Changes the aspect ratio of plots. This might be advanced for some. If so, leave it at default. 
Plotmeinput.axisfontsize    = 26;                                           % Changes the font size of the numbers and labels on the axes. 
Plotmeinput.xtickdecimals   = 1;                                            % Specify how many decimal places you want on your x-axis numbers
Plotmeinput.ytickdecimals   = 2;                                            % Specify how many decimal places you want on your y-axis numbers
Plotmeinput.majorgridoption = 'On';                                         % Major grid on -> 'On' (case sensitive). Major grid off -> Any keyword except 'On'. Use 'Off' to avoid confusion later. 
Plotmeinput.minorgridoption = 'Off';                                        % Minor grid on -> 'On' (case sensitive). Minor grid off -> Any keyword except 'On'. Use 'Off' to avoid confusion later. 
Plotmeinput.legendoption    = 'On';                                         % Want a legend? -> 'On' (case sensitive). Legend off -> Any keyword except 'On'. Use 'Off' to avoid confusion later. 
Plotmeinput.legend          = {'Linear', 'Quadratic',  'Sine',  'Logarithm'};% For n number of data sets, n number of legend entries are required. In this tutorial, n = 4.
Plotmeinput.legendtextsize  = 24;                                           % Specify the legend font size here.                                      
Plotmeinput.legendlocation  = 'best';                                       % Specify legend location. "Best" works fine most times. More more legend location options, consult MATLAB documentation. 
Plotmeinput.textboxoption   = 'Off';                                        % Want to put a text box? -> 'On'. You know by now what to do to leave it at off. Might be advanced. Don't use it if not familiar. 
Plotmeinput.textboxtext     = {'',''};                                      % Specify the text-box text here.
Plotmeinput.textboxdim      = []';                                          % Normalized units. Lower left corner (0, 0) and top right corner (1,1)
Plotmeinput.filename        = 'Example_CoPlot';                             % Specify the file name for the figure. 
Plotmeinput.filelocation    = 'Graphics/Tutorial/Plots';                    % To save the file in the same folder as this MATLAB script use ' ' (needs to have one space between two apostrophes). To save it in a different folder, specify the FilePath. If the folders don't exist, this script will automatically create them.
Plotmeinput.backlocation    = '../../../';                                  % Need to specify the return path to the folder from which this script is running to avoid errors. 
Plotmeinput.nVar            = 4;                                            % Specify the number 'n' here. In this tutorial n = 4, where n is the number of datasets you want to co-plot on a single figure.
Plotme(Plotmeinput);                                                        % Semicolon is required at the end to suppress output from Plotme.m

% The script will save the MATLAB figure in .eps, .emf, .png, .fig, and 
% .pdf file formats in the specified folder location.

%% Credits................................................................

% Created by Suraj Bansal
% PhD Candidate
% Institute for Aerospace Studies
% University of Toronto
% Created: October 10, 2021