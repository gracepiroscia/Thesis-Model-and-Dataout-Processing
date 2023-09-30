% Author: Grace Piroscia
%
% Script to validate triangulation and perform two-way t-tests to determine
% the significance between...
clc;
clear;
%%  Load data
load("data_out_processed.mat")
field_names = {"sc1","sc2","sc3","sc4","sc5","sc6","sc7"};

% Data frames: jth column represents jth scenario. Row i corresponds to
% partipant with id i.
LIKERT_TRUST_RATINGS = [];
MEAN_DTU = []; %(uu)
GAZE_RATIO = [];
TIME_TO_CROSS = []; %(s)
norm_MEAN_DTU = []; 
norm_GAZE_RATIO = [];
norm_TIME_TO_CROSS = []; 
for i = 1:length(field_names)
    LIKERT_TRUST_RATINGS = [LIKERT_TRUST_RATINGS, vertcat(data_processed.(field_names{i}).mean_Likert_rating)];
    MEAN_DTU = [MEAN_DTU,vertcat(data_processed.(field_names{i}).Mean_DTC_uu)];
    GAZE_RATIO = [GAZE_RATIO, vertcat(data_processed.(field_names{i}).Gaze_Ratio)];
    TIME_TO_CROSS = [TIME_TO_CROSS, vertcat(data_processed.(field_names{i}).Cross_Time_s)];
    norm_MEAN_DTU = [norm_MEAN_DTU,vertcat(data_processed.(field_names{i}).norm_Mean_DTC_uu)];
    norm_GAZE_RATIO = [norm_GAZE_RATIO, vertcat(data_processed.(field_names{i}).norm_Gaze_Ratio)];
    norm_TIME_TO_CROSS = [norm_TIME_TO_CROSS, vertcat(data_processed.(field_names{i}).norm_Cross_Time_s)];
end

%% Correlation between subjective and objective measures
% Convert matrix to single col vectors
LIKERT_TRUST_RATINGS_col = LIKERT_TRUST_RATINGS(:);
MEAN_DTU_col = MEAN_DTU(:);
GAZE_RATIO_col = GAZE_RATIO(:);
TIME_TO_CROSS_col = TIME_TO_CROSS(:);
DF = [LIKERT_TRUST_RATINGS_col,MEAN_DTU_col,GAZE_RATIO_col, TIME_TO_CROSS_col];
R = corrcoef(DF);

% plot
figure(1)
set(gcf,'Color', 'w')
subplot(1,3,1)
plot(LIKERT_TRUST_RATINGS_col, MEAN_DTU_col,'o')
ylabel("Mean DTC (uu)")
xlabel("Likert Trust Rating (1-7)")
txt = "R = " + num2str(R(1,2));
text(2,280,txt)
grid on;
subplot(1,3,2)
plot(LIKERT_TRUST_RATINGS_col, GAZE_RATIO_col,'o')
xlabel("Likert Trust Rating (1-7)")
ylabel("Gaze Ratio"); grid on;
txt = "R = " + num2str(R(1,3));
text(2,0.35,txt)
subplot(1,3,3)
plot(LIKERT_TRUST_RATINGS_col, TIME_TO_CROSS_col,'o')
xlabel("Likert Trust Rating (1-7)")
ylabel("Time to Cross (s)")
grid on;
txt = "R = " + num2str(R(1,4));
text(2,16,txt)
A = [LIKERT_TRUST_RATINGS_col,MEAN_DTU_col,GAZE_RATIO_col, TIME_TO_CROSS_col];
R = corrcoef(A);

%% Comparison between eHMI, no eHMI FAV and Human-Driven
% First - lets see how well the eHMI was recognised
eHMI_likert_ratings = vertcat(data_processed.sc2.Likert_eHMI_recognition_rating);
stats_eHMI_likert_rating = [mean(eHMI_likert_ratings), median(eHMI_likert_ratings)];
fprintf("Mean eHMI recognition (Likert): %f\n",stats_eHMI_likert_rating(1));
fprintf("Median eHMI recognition (Likert): %f\n",stats_eHMI_likert_rating(2));

% Whole comparison

% Comparison with just the partipants who had a Likert rating greater than
% the median for recognition of the eHMI


%% Comparison between all non-eHMI FAV scenarios
