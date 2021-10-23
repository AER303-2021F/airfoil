%This code returns scanivalve pressure measurements and Betz pressure
%measurements for 10 velocities and calculates the calibration coefficients
%for the scanivalve
clear all
clc
format long
%%
t_s = 1; %seconds
f_s = 5000; %Hz

prompt  = 'Ensure scanivalve is on channel 0 and the Betz manometer is set to 0, then press enter';
ask     = input(prompt);

%% Communicatng with DAQ...................................................
disp(' *** Initializing the DAQ *** ');

p_scanivalve        = [];
p_Betz              = [];
s                   = daq.createSession('ni');                                  % Creates a session
ch                  = addAnalogInputChannel(s,'Dev1', 1, 'Voltage');            % adds one channel to the session. Standard is differential input
s.DurationInSeconds = t_s;
s.Rate              = f_s;
ch.Range            = [-2.5, 2.5];                                                  % set the channel input range to +/-2.5 Volts

%% Acquiring data...........................................................
Keyword = 'Continue';
key     = 'y';

i = 1;

while strcmp(Keyword,'End') ~= 1
    
    disp([' *** Taking calibration measurement ' num2str(i) ' of 10 *** ']);
    
    if i>1
        
        if strcmp(Keyword, 'End') == 0
            prompt  = 'Change tunnel velocity and press enter'; %time for students to change velocity
            ask     = input(prompt);
            
            disp('Waiting a few seconds for the betz manometer and the wind tunnel to settle')
            pause(15);
            i = i+1;
        end
        
    else
        
        i = i +1;
        
    end

    [data, time]    = s.startForeground;                                       % take some samples
    
    mean_data       = mean(data);
    p_scanivalve(i) = mean_data;
    
    p_Betz(i)       = input('Enter the Betz manometer pressure im mm H20 ');% enter the pressure in mm of H2O
    %p_IM(i)         = input('Enter the inclined manometer pressure im cm h2o ');% enter the pressure in mm of H2O
    %p_Betz(i)       = (10- p_IM(i))*cosd(44)*10;
    key             = input('Press enter, or type "n" to end the calibration','s');
    
    if isempty(key) == 1
        Keyword = 'Continue';
    elseif strcmp(key, 'No') == 1 || strcmp(key, 'no') == 1 || strcmp(key, 'n') == 1
        Keyword = 'End';
    else
        Keyword = 'Continue';
    end
   
    
    
end

disp('Calibration complete');
disp('Calculating calibration coefficients');


%% Calculate calibration coefficients.......................................
coeffs  = polyfit(p_scanivalve(1:end),p_Betz(1:end),1);
x       = linspace(p_scanivalve(1),p_scanivalve(end));
y       = polyval(coeffs,x);

%% Plot calibration curve...................................................
plot(p_scanivalve(1:end),p_Betz(1:end),'xk');
hold on
plot(x,y);
xlabel('p_scanivalve [V]');
ylabel('p_Betz [mm H2O]');

%% Saving calibration data..................................................
disp('Saving calibration data to Calibration_data.mat');
cd Data
save Calibration_data p_scanivalve p_Betz coeffs

delete(s)
clear s