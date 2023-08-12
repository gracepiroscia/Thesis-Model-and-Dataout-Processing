% Author: Grace Piroscia
%
% Script to process the raw .csv data_out files collected during
% participant trials. This file contains the information (in order):
% - Delta Time: time interval between unreal engine event ticks in seconds
% 
% - Player Position: The X,Y,Z position of the HMD in the unreal world in
%   uu units
%
% - Data_out Actor Position: The X,Y,Z position of the center of the
%                            "collision zone". This is constant with time. 
%                            In uu units.
%
% - DTC_uu: The distance between the player and the collision zone with
%           time. In uu units.
%
% - Gaze: The forward vector apporximating gaze as the center of the HMD.
%         In uu units.

clc;
clear;

%% Load data
opts = delimitedTextImportOptions('NumVariables',5);
opts.Delimiter = ",";
opts.VariableNames = {'Delta_secs', 'Player_Loc', 'Collision_Zone_Loc', ...
            'DTC_uu', 'Estimated_Gaze_FV'};
data_out = readtable("csv_files/8_12_11_12.csv", opts);
data_out(1,:) = []; %delete first row with var names

%% Thresholds 
dtc_min_x = 1; %uu
dtc_max_x = 460; %""
gaze_mins_xyz = [200,300,400]; %""
gaze_maxs_xyz = [200,300,400]; % ""

%% Process data
% Convert to struct with relevant data
delta_t_char = data_out.Delta_secs;
play_loc_char = data_out.Player_Loc;
collision_zone_loc_char = data_out.Collision_Zone_Loc;
col_zone_vec = [];
DTC_char = data_out.DTC_uu;
gaze_dir_char = data_out.Estimated_Gaze_FV;

for i = 1:size(data_out,1)

    % Load current data
    delt = str2double(delta_t_char{i});
    playloc = play_loc_char(i);
    dtc = str2double(DTC_char(i));
    gaze_dir = gaze_dir_char(i);

    if i == 1
        s(i).ElapsedTime_s = delt;

        % only need to convert collision zone loc once (static)
        col_loc = collision_zone_loc_char(i);
        xyz_idxs = strfind(col_loc, "=");
        xyz_idxs = xyz_idxs{1,1} + 1;
        space_idxs = strfind(col_loc," ");
        space_idxs = space_idxs{1,1};
        col_loc = col_loc{1,1};
        x = str2double(col_loc( xyz_idxs(1):(space_idxs(1)-1) ));
        y = str2double(col_loc( xyz_idxs(2):(space_idxs(2)-1) ));
        z = str2double(col_loc( xyz_idxs(3):end));
        col_zone_vec = [x,y,z];
        
    else 
        s(i).ElapsedTime_s = s(i-1).ElapsedTime_s + delt;
    end

    % Convert string arrays to vector in the form of [X,Y,Z] positions
    xyz_idxs = strfind(playloc, "=");
    xyz_idxs = xyz_idxs{1,1} + 1;
    space_idxs = strfind(playloc," ");
    space_idxs = space_idxs{1,1};
    playloc = playloc{1,1};
    x = str2double(playloc( xyz_idxs(1):(space_idxs(1)-1) ));
    y = str2double(playloc( xyz_idxs(2):(space_idxs(2)-1) ));
    z = str2double(playloc( xyz_idxs(3):end));
    s(i).PlayerLoc_uu = [x,y,z];

    if (x >= dtc_min_x) && (x <= dtc_max_x)
        s(i).In_DTC_Thresh = 1;
    else
        s(i).In_DTC_Thresh = 0;
    end

    s(i).CollisionZoneLoc_uu = col_zone_vec;

    s(i).DTC_uu = dtc;

    xyz_idxs = strfind(gaze_dir, "=");
    xyz_idxs = xyz_idxs{1,1} + 1;
    space_idxs = strfind(gaze_dir," ");
    space_idxs = space_idxs{1,1};
    gaze_dir = gaze_dir{1,1};
    x = str2double(gaze_dir( xyz_idxs(1):(space_idxs(1)-1) ));
    y = str2double(gaze_dir( xyz_idxs(2):(space_idxs(2)-1) ));
    z = str2double(gaze_dir( xyz_idxs(3):end));
    s(i).EstimatedGazeVec_uu = [x,y,z];

    if (x>= gaze_mins_xyz(1) && x<= gaze_maxs_xyz(1)) && (y>= gaze_mins_xyz(2)...
            && y<= gaze_maxs_xyz(2)) && (z>= gaze_mins_xyz(3) && z<= gaze_maxs_xyz(3))
        s(i).Is_Looking_at_Vehicle = 1;
    else 
        s(i).Is_Looking_at_Vehicle = 0;
    end
   
end

%% Plot data
% DTC and trajectory with time
bool_InColZone = vertcat(s.In_DTC_Thresh);
df = [vertcat(s.ElapsedTime_s), vertcat(s.DTC_uu), vertcat(s.PlayerLoc_uu)];
df = df(bool_InColZone == 1,:); %keep the relevant rows (to start averaging collision)
avg_DTC = mean(df(:,2));
fprintf("Average DTC = %f\n", avg_DTC)

figure(1)
set(gcf,'Color','w');
subplot(1,2,1)
plot(df(:,1), df(:,2)); grid on;
xlabel('Time (s)')
ylabel('Distance to Collision Zone (uu)')
subplot(1,2,2)
col_coord = s(1).CollisionZoneLoc_uu; 
plot(col_coord(1), col_coord(2), 'ko', 'MarkerSize', 5)
hold on
plot(df(:,3), df(:,4)); grid on;
xlabel("X Position (uu)")
ylabel("Y Position (uu)")

% Boolean Gaze 
bool_Gaze = vertcat(s.Is_Looking_at_Vehicle);
time = vertcat(s.ElapsedTime_s);
figure(2)
set(gcf,'Color','w');
plot(time, bool_Gaze, '*r');
xlabel("Time (s)")
ylabel("Bool Looking at Vehicle (1=yes, 0=no)")

