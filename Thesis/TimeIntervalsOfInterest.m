% Author: Grace Piroscia
%
% Script to go through each of the relevant structs produced by
% ProcessBag.m to determine the time intervals of interest. These will be
% where the buggy has come to a pedestrian crossing with the behaviour
% broken down into to different interactions:
%   1. Complete stop
%   2. Crawl
%
% These times of interest are able to be determined by viewing the
% corresponding buggy camera footage and correlating certain events to the
% motion of the buggy, as plotted by this script.
%
% Looking for non-slope crossings either with multiple or single
% pedestrian.

clc;
clear;

%% Load relevant bag data
FileName = '2018-04-19-15-05-55_Dataset_Year.mat';
load(FileName);

%% Plotting 
AX = vertcat(s.imu.LinearAcc_X);
AZ = vertcat(s.imu.LinearAcc_Z);
time = vertcat(s.imu.Time_since_start_secs);

% Acceleration profiles
figure(1)
subplot(2,1,1)
set(gcf, 'color', 'w');
plot(time,AX, '.k'); grid on;
title('Linear Acceleration - X', 'FontSize',16)
ylabel('Acceleration (m/s^2)', 'FontSize',15)
xlabel('Elapsed Time (s)', 'FontSize',15)
subplot(2,1,2)
plot(time,AZ, '.k'); grid on;
title('Linear Acceleration - Z', 'FontSize',16)

% Velocity Profiles
VX = vertcat(s.odomRear.LinearVel_X);
time = vertcat(s.odomRear.Time_since_start_secs);
figure(2)
set(gcf, 'color', 'w');
plot(time,VX, '.k'); grid on;
title('Linear Velocity - X', 'FontSize',16)
ylabel('Velocity (m/s)', 'FontSize',15)
xlabel('Elapsed Time (s)', 'FontSize',15)


%% Store relevant time intervals of interest. 
%(done visually, with buggy footage and plots to correlate)
% For single pedestrian interactions
CROSSING_STOP_TIME = [71.53, 347.83;82.03,349.83];% Matrix where each column corresponds to the time interval (in time since start)
                     % that the event occured. Top row is the start time,
                     % second row is the end time. (for complete stop
                     % events)
CROSSING_CRAWL_TIME = [57.03,1191.73;64.63,1199.03];%"" but for crawling events

is_Single_Pedestrian_stop = [0,1];% bool array corresponding to each time event, indicating whether the interaction
                            % was between a single pedestrian and the buggy
                            % (1) or multiple pedestrians (0)
                            
is_Single_Pedestrian_crawl = [0,0];

%% Store in struct
s.odomRear(1).Crossing_Complete_Stop_event_timeSinceStart = CROSSING_STOP_TIME;
s.odomRear(1).Crossing_Crawl_event_timeSinceStart = CROSSING_CRAWL_TIME;
s.odomRear(1).IsSinglePedestrianStop = is_Single_Pedestrian_stop;
s.odomRear(1).IsSinglePedestrianStopCrawl = is_Single_Pedestrian_crawl;
save(FileName, 's')
