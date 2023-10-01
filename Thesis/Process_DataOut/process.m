% Author: Grace Piroscia
%
% Script to process the data_out struct, calculating the average dtc and 
% the time taken to cross in the zone of interest as well as gaze ratio. 
% Definitions:
    % - Average Distance-To-Collision: The average distance between the
    %   collision zone location (center of the boundary in the lane of the
    %   vehicle) and the particpant (uu). This is only calculated over the lane
    %   where the vehicle is (zone of interest)
    %
    % - Crossing time: duration in seconds spent by the partipant to cross
    %   the zone of interest.
    %
    % - Gaze Ratio: The time spent looking at the vehicle divded by the entire
    %   time to complete simulation
%
% Note: the normalised version of this data is also calculated (based on
% individual participant maximum and minimum)


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

        % Calculate the time taken to cross zone of interest
        time_in_zone = time(Is_in_DTC_Thresh == 1);
        time_to_cross = time_in_zone(end) - time_in_zone(1);

        % Save to new, processed, struct
        data_processed.(sc_labels{scenario_n})(participant_n).Mean_DTC_uu = mean_dtc;
        data_processed.(sc_labels{scenario_n})(participant_n).Gaze_Ratio = gaze_ratio;
        data_processed.(sc_labels{scenario_n})(participant_n).Cross_Time_s = time_to_cross;

        if gaze_ratio == 0
            disp([participant_n,sc_labels{scenario_n} ]);
        end
    end
end

%% Calculate normalised variables
for participant_n = 1:30
    % Across all scenarios (i.e. 7 values for each)
    mean_dtc = [];
    gaze_ratio = [];
    cross_time =[];
    likert_rating = [];
    for i = 1:7
        mean_dtc = [mean_dtc, data_processed.(sc_labels{i})(participant_n).Mean_DTC_uu];
        gaze_ratio = [gaze_ratio, data_processed.(sc_labels{i})(participant_n).Gaze_Ratio];
        cross_time =[cross_time, data_processed.(sc_labels{i})(participant_n).Cross_Time_s];
        likert_rating = [likert_rating, data_processed.(sc_labels{i})(participant_n).mean_Likert_rating];
    end

    % normalise
    mean_dtc_norm = (mean_dtc - min(mean_dtc))./(max(mean_dtc) - min(mean_dtc));
    gaze_ratio_norm = (gaze_ratio - min(gaze_ratio))./(max(gaze_ratio) - min(gaze_ratio));
    cross_time_norm =(cross_time - min(cross_time))./(max(cross_time) - min(cross_time));
    likert_rating_norm = (likert_rating - min(likert_rating))./(max(likert_rating) - min(likert_rating));

    % Store
    for i = 1:7
        data_processed.(sc_labels{i})(participant_n).norm_Mean_DTC_uu = mean_dtc_norm(i);
        data_processed.(sc_labels{i})(participant_n).norm_Gaze_Ratio = gaze_ratio_norm(i);
        data_processed.(sc_labels{i})(participant_n).norm_Cross_Time_s = cross_time_norm(i);
        data_processed.(sc_labels{i})(participant_n).norm_Mean_Likert_rating = likert_rating_norm(i);
    end

end

%% Save Struct
save("data_out_processed.mat", 'data_processed')
