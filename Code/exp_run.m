%This code returns the surface pressure measurements around the airfoil and
%the wake measurements

clear all
clc
%%
%Enter the sampling time and sampling frequency determined in the first
%part of this lab along with the current angle of attack

t_s = 1;                                                                    % seconds
f_s = 5000;                                                                 % Hz
AoA = 3;                                                                    % Angle of attack

%%
prompt  = 'Make sure the angle of attack in this code has been updated and the scanivalve is on channel 1, then press enter';
ask     = input(prompt);

x = [0.00, 0.03, 0.06, 0.10, 0.15, 0.20, 0.30, 0.40, 00.55, 00.70, 00.85, 01.00, 00.90, 00.60, 00.40, 00.30, 00.20, 00.10, 00.05];
y = [0.00, 1.67, 3.33, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.67, 18.33, 20];
y2 = y + 0.5;

disp(' *** Initializing the DAQ *** ');

s                       = daq.createSession('ni');                                          % Creates a session
ch                      = addAnalogInputChannel(s,'Dev1', 1, 'Voltage');                    % adds one channel to the session. Standard is differential input
s.DurationInSeconds     = t_s;                                                              % Sampling time (time for which the data is acquired)
s.Rate                  = f_s;                                                              % Sampling frequency (the rate at which data is acquired, measured in Hz)
ch.Range                = [-2.5, 2.5];      


[data, time]            = s.startForeground;                                       % take some samples
mean_data               = mean(data);
p_airfoil_tare_V        = mean_data;


% set the channel input range to +/-2.5 Volts. This means, the DAQ channels are capable
                                                                                            % of measuing up to +/- 2.5 Volts from the device which is connected to it.

% %% Take Tare Measurements..................................................
% 
% for i = 1 : 19
%     
%     disp([' *** Taking zero-offset surface pressure measurement ' num2str(i) ' of 19 *** ']);
%     
%     [data, time]            = s.startForeground;                                       % take some samples
%     mean_data               = mean(data);
%     p_airfoil_tare_V(i)     = mean_data;
% 
%     if i < 19
%         prompt  = ['Change to channel ',num2str(i+1), ' and press enter'];              %time for students to change scanivalve channel
%         ask     = input(prompt);
%     end
%     
% end
% 
% fprintf('\n')
% prompt  = ['Change to channel ',num2str(i+1), ' and press enter'];
% ask = input(prompt);
% fprintf('\n')
% 
% for i = 1   :   17
% 
%     disp([' *** Taking zero-offset wake measurements ' num2str(i) ' of 17 *** ']);
%     
%     [data, time]        = s.startForeground; 
%     mean_data           = mean(data);
%     p_rake_tare_V(i)    = mean_data;
%     e_std_rake1(i)      = std(data);
%     
%     if i < 17
%         prompt  = ['Change to channel ', num2str(i+20), ' and press enter']; %time for students to change scanivalve channel
%         ask     = input(prompt);
%     end
%     
% end
% 
% fprintf('\nZero-offset measurements completed\n\n')

Re      = input('Enter the desired Reynolds number: ');
V_mps   = Re * 1.4607e-5 / 0.1;

ask = input(['Set the wind tunnel speed to ', num2str(V_mps), 'm/s, and then press enter']);
fprintf('Settling...\n')
pause(5);
ask = input('Set the Scanivalve channel to 1 and press enter to begin taking data ');
fprintf('\n');

for i = 1 : 19
    disp([' *** Taking surface pressure measurement ' num2str(i) ' of 19 *** ']);
    
    [data, time]    = s.startForeground;                                       % Take some samples at the sampling rate and for the duration of time (sampling time) set before. 
    p_airfoil(i)    = mean(data);                                              % Take mean of those samples

    if i < 19
        prompt  = ['Change to channel ',num2str(i+1), ' and press enter'];     % Time for students to change scanivalve channel. How many of you hate the scanivalve button like I do?
        ask     = input(prompt);
    end
end
fprintf('\nSurface pressure measurements complete\n\n');

prompt  = 'Change to channel 20 and press enter';                               %time for students to change scanivalve channel
ask = input(prompt);

for i = 1   :   17
    disp([' *** Taking wake measurement ' num2str(i) ' of 17 *** ']);
    
    [data, time]    = s.startForeground; 
    p_rake1(i)       = mean(data);
    
    if i < 17
        prompt  = ['Change to channel ', num2str(i+20), ' and press enter'];  %time for students to change scanivalve channel
        ask     = input(prompt);
    end
    
end

fprintf('Wake measurements complete\n');

prompt = 'Traverse wake rake 0.5cm and set the scanivalve to channel 20, then press enter'; %time for students to change scanivalve channel
ask = input(prompt);

k = 0;
for i = 1:17
    disp([' *** Taking wake measurement ' num2str(i) ' of 17 *** ']);
    [data, time] = s.startForeground; % take some samples
    p_rake2(i) = mean(data);
    
    if i < 17
        prompt  = ['Change to channel ', num2str(i+20), ' and press enter']; %time for students to change scanivalve channel
        ask     = input(prompt);
    end
    
end
disp('Run complete');

cd Data
save(['Experimental_data_' mat2str(AoA) '.mat']);
cd ../

delete(s)
clear s
