% Author: Grace Piroscia
%
% Script to process the data_out struct, calculating the average dtc in the
% zone of interest as well as gaze ratio. Definitions:
%
    % - Average Distance-To-Collision: The average distance between the
    %   collision zone location (center of the boundary in the lane of the
    %   vehicle) and the particpant (uu). This is only calculated over the lane
    %   where the vehicle is.
    %
    % - Gaze Ratio: The time spent looking at the vehicle divded by the entire
    %   time to complete simulation

clc
clear
 %% Load data
 load("data_out.mat")
 data_processed = data; 

 %% Calculate avg DTC and gaze ratio
 sc_labels = {"sc1","sc2","sc3","sc4","sc5","sc6","sc7"};

for participant_n = 1:30
    for scenario_n = 1:7
        % Grab relevant variables
        DTC = data.(sc_labels{scenario_n})(participant_n).DTC_uu;
        Is_in_DTC_Thresh = data.(sc_labels{scenario_n})(participant_n).In_DTC_Thresh;
        time = data.(sc_labels{scenario_n})(participant_n).Elapsed_Time_s;
        Is_looking_at_Vehicle =data.(sc_labels{scenario_n})(participant_n).Is_Looking_at_Vehicle;
        
        % Calculate average DTC in zone of interest (thresholds defined in
        % csv_import.m)
        dtc_cut = DTC(Is_in_DTC_Thresh == 1);
        mean_dtc = mean(dtc_cut);

        % Calculate Gaze Ratio
        total_counts_looking = sum(Is_looking_at_Vehicle == 1); %discretising
        total_time_counts = length(time);
        gaze_ratio = total_counts_looking/total_time_counts;

        % Save to new, processed, struct
        data_processed.(sc_labels{scenario_n})(participant_n).Mean_DTC_uu = mean_dtc;
        data_processed.(sc_labels{scenario_n})(participant_n).Gaze_Ratio = gaze_ratio;

        if gaze_ratio == 0
            disp([participant_n,sc_labels{scenario_n} ]);
        end
    end
end

%% Save Struct
save("data_out_processed.mat", 'data_processed')
