% Author: Grace Piroscia
%
% Script to validate triangulation and perform two-sample t-tests to determine
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

% Whole comparison: sc2 = eHMI, sc3 = human-driven, sc4 = no eHMI FAV
mean_Likert = [mean(LIKERT_TRUST_RATINGS(:,2)),mean(LIKERT_TRUST_RATINGS(:,3)),mean(LIKERT_TRUST_RATINGS(:,4))];
mean_Mean_DTC = [mean(MEAN_DTU(:,2)),mean(MEAN_DTU(:,3)),mean(MEAN_DTU(:,4))];
mean_Cross_time = [mean(TIME_TO_CROSS(:,2)),mean(TIME_TO_CROSS(:,3)),mean(TIME_TO_CROSS(:,4))];
mean_Gaze_ratio = [mean(GAZE_RATIO(:,2)),mean(GAZE_RATIO(:,3)),mean(GAZE_RATIO(:,4))];
SEM_Likert = [SEM(LIKERT_TRUST_RATINGS(:,2)),SEM(LIKERT_TRUST_RATINGS(:,3)),SEM(LIKERT_TRUST_RATINGS(:,4))];
SEM_Mean_DTC = [SEM(MEAN_DTU(:,2)),SEM(MEAN_DTU(:,3)),SEM(MEAN_DTU(:,4))];
SEM_Cross_time = [SEM(TIME_TO_CROSS(:,2)),SEM(TIME_TO_CROSS(:,3)),SEM(TIME_TO_CROSS(:,4))];
SEM_Gaze_ratio = [SEM(GAZE_RATIO(:,2)),SEM(GAZE_RATIO(:,3)),SEM(GAZE_RATIO(:,4))];


% T-test
for i = 1:3
    for j = 1:3
        [h_likert(i,j), p_likert(i,j)]  = ttest2(LIKERT_TRUST_RATINGS(:,(i+1)), LIKERT_TRUST_RATINGS(:,(j+1)));
        [h_meanDTC(i,j), p_meanDTC(i,j)]  = ttest2(MEAN_DTU(:,(i+1)), MEAN_DTU(:,(j+1)));
        [h_crossTime(i,j), p_crossTime(i,j)]  = ttest2(TIME_TO_CROSS(:,i+1), TIME_TO_CROSS(:,j+1));
        [h_gazeRatio(i,j), p_gazeRatio(i,j)]  = ttest2(GAZE_RATIO(:,i+1), GAZE_RATIO(:,j+1));
    end
end
% ^ No sig diff. Closest is for Likert between human driven and no eHMI FAV at p = 0.052     

% Plot
x = 1:3;
x_labels = {"eHMI FAV", "Human-Driven", "no eHMI FAV"} ;
figure('Name','Differences in eHMI - all')
set(gcf,'Color', 'w'); set(gcf, 'Position',  [100, 100, 1100, 700]); subplot(2,2,1)
bar(x,mean_Likert); hold on;
e = errorbar(x,mean_Likert, SEM_Likert); grid on; e.Color = 'black'; e.LineWidth = 1;
set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
ylabel("Likert Subjective-Trust Score (1-7)", 'FontSize',14)
subplot(2,2,2)
bar(x,mean_Mean_DTC); hold on;
e= errorbar(x,mean_Mean_DTC, SEM_Mean_DTC); grid on;e.Color = 'black'; e.LineWidth = 1;
ylabel("Mean DTC (uu)", 'FontSize',14); ylim([220,260]);
set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
subplot(2,2,3)
bar(x,mean_Cross_time); hold on;
e= errorbar(x,mean_Cross_time, SEM_Cross_time); grid on;e.Color = 'black'; e.LineWidth = 1;
ylabel("Time to Cross DTC-Zone (s)", 'FontSize',14)
set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
subplot(2,2,4)
bar(x,mean_Gaze_ratio); hold on;
e=errorbar(x,mean_Gaze_ratio, SEM_Gaze_ratio); grid on;e.Color = 'black'; e.LineWidth = 1;
set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
ylabel("Gaze Ratio", 'FontSize',14)



% Comparison with just the partipants who had a Likert rating greater than
% the median for recognition of the eHMI


%% Comparison between all non-eHMI FAV scenarios
