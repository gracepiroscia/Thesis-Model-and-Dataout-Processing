% Author: Grace Piroscia
%
% Script to process the raw .csv files collected during
% participant trials for both the objective (unreal engine measurements) 
% and subjective (Likert-Scale questionnaire) data, collating them into a
% single struct.
% 
% Each scenario file per participant for the OBJECTIVE data contains the 
% information (in order):
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
%
% Each scenario file per participant for the SUBJECTIVE data contains the
% information (in order):
    %  - The Likert score (1-7) per subjective-trust question (4), per
    %    scenario, identical to all scenarios
    %
    %  - The average Likert score (1-7) across all four subjective-trust
    %    questions
    % 
    %  - (sc1,2,4-7) Boolean indicator representing if particpants 
    %    regonised there WAS NOT a driver (0 = no, 1 = yes)
    %
    %  - (sc3) Boolean indicator representing if particpants regonised there 
    %     WAS a driver (0 = no, 1 = yes)
    %
    %  - (sc2 only) Likert score (1-7) representing how well the
    %     particpant noticed the eHMI design.

clc;
clear;

%% Plot data option
% Will plot per participant per scenario the DTC with time and gaze bool
% with time. Will wait for user to press any key to continue to next
% scenario.
plot_data = false;


%% Thresholds 
dtc_min_x = 10; %uu
dtc_max_x = 600; %""
gaze_mins_xyz = [15,-100,-60]; %""
gaze_maxs_xyz = [70,-50,20];%[50,-80,0]; % ""

%% 1. Load OBJECTIVE, unreal engine data
opts = delimitedTextImportOptions('NumVariables',5);
opts.Delimiter = ",";
opts.VariableNames = {'Delta_secs', 'Player_Loc', 'Collision_Zone_Loc', ...
            'DTC_uu', 'Estimated_Gaze_FV'};
sc_labels = {"sc1","sc2","sc3","sc4","sc5","sc6","sc7"};

for participant_n = 1:30
dir_name = "UE_data/p" + num2str(participant_n);
FileNames = {dir(dir_name).name};
    for scenario_n = 1:7
        str_of_interest = sc_labels{scenario_n} + ".csv";
        file_name = "";

        % Find the csv file relevant to current scenario
        for i = 1:length(FileNames)
            str_name = convertCharsToStrings(FileNames{1,i});
            k = strfind(str_name, str_of_interest);
            if ~isempty(k)
                file_name = str_name;
            end
        end

        %Load csv file
        full_fileName = dir_name + "/"+ file_name;
        data_out = readtable(full_fileName, opts);
        data_out(1,:) = []; %delete first row with var names6
        data.(sc_labels{scenario_n})(participant_n).p_id = participant_n;

        %% Process data
        % Convert to struct with relevant data
        s = []; % Temporary struct to load in data
        delta_t_char = data_out.Delta_secs;
        play_loc_char = data_out.Player_Loc;
        collision_zone_loc_char = data_out.Collision_Zone_Loc;
        col_zone_vec = [];
        DTC_char = data_out.DTC_uu;
        gaze_dir_char = data_out.Estimated_Gaze_FV;

        for i = 1:size(data_out,1) %for each row in csv file

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
                if (col_zone_vec(1) ~= 200.00) || (col_zone_vec(2) ~= -480.00) || (col_zone_vec(3) ~= 60.00)
                    disp(full_fileName);
                end


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

        % Store to data struct
        data.(sc_labels{scenario_n})(participant_n).Elapsed_Time_s = vertcat(s.ElapsedTime_s);
        data.(sc_labels{scenario_n})(participant_n).Play_Loc_xyz_uu = vertcat(s.PlayerLoc_uu);
        data.(sc_labels{scenario_n})(participant_n).DTC_uu = vertcat(s.DTC_uu);
        data.(sc_labels{scenario_n})(participant_n).In_DTC_Thresh = vertcat(s.In_DTC_Thresh);
        data.(sc_labels{scenario_n})(participant_n).EstimatedGazeVec_xyz_uu = vertcat(s.EstimatedGazeVec_uu);
        data.(sc_labels{scenario_n})(participant_n).Is_Looking_at_Vehicle = vertcat(s.Is_Looking_at_Vehicle);

        % Option to plot
        if plot_data 
            % DTC and trajectory with time
            bool_InColZone = vertcat(s.In_DTC_Thresh);
            df = [vertcat(s.ElapsedTime_s), vertcat(s.DTC_uu), vertcat(s.PlayerLoc_uu)];
            df = df(bool_InColZone == 1,:); %keep the relevant rows (to start averaging collision)
            avg_DTC = mean(df(:,2));
            %fprintf("Average DTC = %f\n", avg_DTC)

            f1 = figure(1);
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
            xyz_gaze = vertcat(s.EstimatedGazeVec_uu);
            f2 = figure(2);
            set(gcf,'Color','w');
            subplot(1,2,1)
            plot(time, bool_Gaze, '*r');
            xlabel("Time (s)")
            ylabel("Bool Looking at Vehicle (1=yes, 0=no)")
            subplot(1,2,2)
            hold on
            plot(time, xyz_gaze(:,1))
            plot(time, xyz_gaze(:,2))
            plot(time, xyz_gaze(:,3))
            legend("x","y","z")

            w = waitforbuttonpress;
            clf(f1);
            clf(f2);
        end 

    end
end

%% 2. Load SUBJECTIVE, likert questionnaire data
SheetNames = {"sc1","sc2","sc3","sc4","sc5","sc6","sc7"};

for participant_n = 1:30
    for scenario_n = 1:7
        CurrentScenarioSheet = readtable("latin_square_scenario.xlsx", 'Sheet',SheetNames{scenario_n});
        CurrentScenarioSheet(1:2,:) = []; %delete first + second row 
        avg_Likert_rating = CurrentScenarioSheet.Var6(participant_n);
        data.(SheetNames{scenario_n})(participant_n).mean_Likert_rating = avg_Likert_rating;

        if SheetNames{scenario_n} ~= "sc3"
            data.(SheetNames{scenario_n})(participant_n).noticed_no_driver = ...
                CurrentScenarioSheet.Var7(participant_n);
        else
            data.(SheetNames{scenario_n})(participant_n).noticed_driver = ...
                CurrentScenarioSheet.Var7(participant_n);
        end

        if SheetNames{scenario_n} == "sc2"
            data.(SheetNames{scenario_n})(participant_n).Likert_eHMI_recognition_rating = ...
                CurrentScenarioSheet.Var8(participant_n);        
        end
    end
end

%% Save struct
save("data_out.mat", 'data');





