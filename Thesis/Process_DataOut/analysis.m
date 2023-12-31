% Author: Grace Piroscia
%
% Script to validate triangulation and perform two-sample t-tests to determine
% the significance between...
clc;
clear;
%%  Load data
load("Process_DataOut/data_out_processed.mat")
field_names = {"sc1","sc2","sc3","sc4","sc5","sc6","sc7"};

% Data frames: jth column represents jth scenario. Row i corresponds to
% partipant with id i.
LIKERT_TRUST_RATINGS = [];
MEAN_DTU = []; %(uu)
GAZE_RATIO = [];
TIME_TO_CROSS = []; %(s)
norm_LIKERT_TRUST_RATINGS = [];
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
    norm_LIKERT_TRUST_RATINGS = [norm_LIKERT_TRUST_RATINGS, vertcat(data_processed.(field_names{i}).norm_Mean_Likert_rating)];
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
ylabel("Mean DTC (uu)", 'FontSize',14)
xlabel("Likert Trust Rating (1-7)", 'FontSize',14)
txt = "R = " + num2str(R(1,2));
t = text(2,280,txt);
t.FontSize = 14;
grid on;
subplot(1,3,2)
plot(LIKERT_TRUST_RATINGS_col, GAZE_RATIO_col,'o')
xlabel("Likert Trust Rating (1-7)", 'FontSize',14)
ylabel("Gaze Ratio", 'FontSize',14); grid on;
title("Comparison between Pedestrian Metrics", 'FontSize',16)
txt = "R = " + num2str(R(1,3));
t = text(2,0.35,txt);
t.FontSize = 14;
subplot(1,3,3)
plot(LIKERT_TRUST_RATINGS_col, TIME_TO_CROSS_col,'o')
xlabel("Likert Trust Rating (1-7)", 'FontSize',14)
ylabel("Time to Cross (s)", 'FontSize',14)
grid on;
txt = "R = " + num2str(R(1,4));
t = text(2,16,txt);
t.FontSize = 14;
A = [LIKERT_TRUST_RATINGS_col,MEAN_DTU_col,GAZE_RATIO_col, TIME_TO_CROSS_col];
R = corrcoef(A);

%% Comparison between eHMI, no eHMI FAV and Human-Driven - NON-normalised
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
% the median ( > 5 i.e 6 or more) for recognition of the eHMI
% Note: we have N = 14
id_recogise_eHMI = find(eHMI_likert_ratings>5);
%filter original arrays
mean_Likert_filt = [mean(LIKERT_TRUST_RATINGS(id_recogise_eHMI,2)),mean(LIKERT_TRUST_RATINGS(id_recogise_eHMI,3)),mean(LIKERT_TRUST_RATINGS(id_recogise_eHMI,4))];
mean_Mean_DTC_filt = [mean(MEAN_DTU(id_recogise_eHMI,2)),mean(MEAN_DTU(id_recogise_eHMI,3)),mean(MEAN_DTU(id_recogise_eHMI,4))];
mean_Cross_time_filt = [mean(TIME_TO_CROSS(id_recogise_eHMI,2)),mean(TIME_TO_CROSS(id_recogise_eHMI,3)),mean(TIME_TO_CROSS(id_recogise_eHMI,4))];
mean_Gaze_ratio_filt = [mean(GAZE_RATIO(id_recogise_eHMI,2)),mean(GAZE_RATIO(id_recogise_eHMI,3)),mean(GAZE_RATIO(id_recogise_eHMI,4))];
SEM_Likert_filt = [SEM(LIKERT_TRUST_RATINGS(id_recogise_eHMI,2)),SEM(LIKERT_TRUST_RATINGS(id_recogise_eHMI,3)),SEM(LIKERT_TRUST_RATINGS(id_recogise_eHMI,4))];
SEM_Mean_DTC_filt = [SEM(MEAN_DTU(id_recogise_eHMI,2)),SEM(MEAN_DTU(id_recogise_eHMI,3)),SEM(MEAN_DTU(id_recogise_eHMI,4))];
SEM_Cross_time_filt = [SEM(TIME_TO_CROSS(id_recogise_eHMI,2)),SEM(TIME_TO_CROSS(id_recogise_eHMI,3)),SEM(TIME_TO_CROSS(id_recogise_eHMI,4))];
SEM_Gaze_ratio_filt = [SEM(GAZE_RATIO(id_recogise_eHMI,2)),SEM(GAZE_RATIO(id_recogise_eHMI,3)),SEM(GAZE_RATIO(id_recogise_eHMI,4))];

% T-test
for i = 1:3
    for j = 1:3
        [h_likert_f(i,j), p_likert_f(i,j)]  = ttest2(LIKERT_TRUST_RATINGS(id_recogise_eHMI,(i+1)), LIKERT_TRUST_RATINGS(id_recogise_eHMI,(j+1)));
        [h_meanDTC_f(i,j), p_meanDTC_f(i,j)]  = ttest2(MEAN_DTU(id_recogise_eHMI,(i+1)), MEAN_DTU(id_recogise_eHMI,(j+1)));
        [h_crossTime_f(i,j), p_crossTime_f(i,j)]  = ttest2(TIME_TO_CROSS(id_recogise_eHMI,i+1), TIME_TO_CROSS(id_recogise_eHMI,j+1));
        [h_gazeRatio_f(i,j), p_gazeRatio_f(i,j)]  = ttest2(GAZE_RATIO(id_recogise_eHMI,i+1), GAZE_RATIO(id_recogise_eHMI,j+1));

    end
end
% ^ Sig diff between Likert score for no eHMI FAV and human-driven 
fprintf("Statistical Significance (Likert): Filtered by recognition, Human driven vs no eHMI FAV, p = %f\n", p_likert_f(2,3))

% Plot
x = 1:3;
x_labels = {"eHMI FAV", "Human-Driven", "no eHMI FAV"} ;
figure('Name','Differences in eHMI - Likert Recognition > 5')
set(gcf,'Color', 'w'); set(gcf, 'Position',  [100, 100, 1100, 700]); subplot(2,2,1)
bar(x,mean_Likert_filt); hold on;
e = errorbar(x,mean_Likert_filt, SEM_Likert_filt); grid on; e.Color = 'black'; e.LineWidth = 2;
set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
ylabel("Likert Subjective-Trust Score (1-7)", 'FontSize',14)
title("Trust in Vehicle", 'FontSize',16)
ylim([0,7])
subplot(2,2,2)
bar(x,mean_Mean_DTC_filt); hold on;
e= errorbar(x,mean_Mean_DTC_filt, SEM_Mean_DTC_filt); grid on;e.Color = 'black'; e.LineWidth = 1;
ylabel("Mean DTC (uu)", 'FontSize',14); ylim([220,260]);
set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
subplot(2,2,3)
bar(x,mean_Cross_time_filt); hold on;
e= errorbar(x,mean_Cross_time_filt, SEM_Cross_time_filt); grid on;e.Color = 'black'; e.LineWidth = 1;
ylabel("Time to Cross DTC-Zone (s)", 'FontSize',14)
set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
subplot(2,2,4)
bar(x,mean_Gaze_ratio_filt); hold on;
e=errorbar(x,mean_Gaze_ratio_filt, SEM_Gaze_ratio_filt); grid on;e.Color = 'black'; e.LineWidth = 1;
set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
ylabel("Gaze Ratio", 'FontSize',14)

%% Comparison between eHMI, no eHMI FAV and Human-Driven - NORMALISED 
%% MAKES NO DIFF TO SIGNIFICANCE
% % Whole comparison: sc2 = eHMI, sc3 = human-driven, sc4 = no eHMI FAV
% mean_Likert = [mean(norm_LIKERT_TRUST_RATINGS(:,2)),mean(norm_LIKERT_TRUST_RATINGS(:,3)),mean(norm_LIKERT_TRUST_RATINGS(:,4))];
% mean_Mean_DTC = [mean(norm_MEAN_DTU(:,2)),mean(norm_MEAN_DTU(:,3)),mean(norm_MEAN_DTU(:,4))];
% mean_Cross_time = [mean(norm_TIME_TO_CROSS(:,2)),mean(norm_TIME_TO_CROSS(:,3)),mean(norm_TIME_TO_CROSS(:,4))];
% mean_Gaze_ratio = [mean(norm_GAZE_RATIO(:,2)),mean(norm_GAZE_RATIO(:,3)),mean(norm_GAZE_RATIO(:,4))];
% SEM_Likert = [SEM(norm_LIKERT_TRUST_RATINGS(:,2)),SEM(norm_LIKERT_TRUST_RATINGS(:,3)),SEM(norm_LIKERT_TRUST_RATINGS(:,4))];
% SEM_Mean_DTC = [SEM(norm_MEAN_DTU(:,2)),SEM(norm_MEAN_DTU(:,3)),SEM(norm_MEAN_DTU(:,4))];
% SEM_Cross_time = [SEM(norm_TIME_TO_CROSS(:,2)),SEM(norm_TIME_TO_CROSS(:,3)),SEM(norm_TIME_TO_CROSS(:,4))];
% SEM_Gaze_ratio = [SEM(norm_GAZE_RATIO(:,2)),SEM(norm_GAZE_RATIO(:,3)),SEM(norm_GAZE_RATIO(:,4))];
% 
% % T-test
% for i = 1:3
%     for j = 1:3
%         [h_likert_norm(i,j), p_likert_norm(i,j)]  = ttest2(norm_LIKERT_TRUST_RATINGS(:,(i+1)), norm_LIKERT_TRUST_RATINGS(:,(j+1)));
%         [h_meanDTC_norm(i,j), p_meanDTC_norm(i,j)]  = ttest2(norm_MEAN_DTU(:,(i+1)), norm_MEAN_DTU(:,(j+1)));
%         [h_crossTime_norm(i,j), p_crossTime_norm(i,j)]  = ttest2(norm_TIME_TO_CROSS(:,i+1), norm_TIME_TO_CROSS(:,j+1));
%         [h_gazeRatio_norm(i,j), p_gazeRatio_norm(i,j)]  = ttest2(norm_GAZE_RATIO(:,i+1), norm_GAZE_RATIO(:,j+1));
%     end
% end
% 
% % Plot
% x = 1:3;
% x_labels = {"eHMI FAV", "Human-Driven", "no eHMI FAV"} ;
% figure('Name','Differences in eHMI - all, NORMALISED')
% set(gcf,'Color', 'w'); set(gcf, 'Position',  [100, 100, 1100, 700]); subplot(2,2,1)
% bar(x,mean_Likert); hold on;
% e = errorbar(x,mean_Likert, SEM_Likert); grid on; e.Color = 'black'; e.LineWidth = 1;
% set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
% ylabel("Likert Subjective-Trust Score (1-7)", 'FontSize',14)
% subplot(2,2,2)
% bar(x,mean_Mean_DTC); hold on;
% e= errorbar(x,mean_Mean_DTC, SEM_Mean_DTC); grid on;e.Color = 'black'; e.LineWidth = 1;
% ylabel("Mean DTC (uu)", 'FontSize',14);
% set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
% subplot(2,2,3)
% bar(x,mean_Cross_time); hold on;
% e= errorbar(x,mean_Cross_time, SEM_Cross_time); grid on;e.Color = 'black'; e.LineWidth = 1;
% ylabel("Time to Cross DTC-Zone (s)", 'FontSize',14)
% set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
% subplot(2,2,4)
% bar(x,mean_Gaze_ratio); hold on;
% e=errorbar(x,mean_Gaze_ratio, SEM_Gaze_ratio); grid on;e.Color = 'black'; e.LineWidth = 1;
% set(gca,'xticklabel',x_labels, 'fontsize',14); xtickangle(gca,30);
% ylabel("Gaze Ratio", 'FontSize',14)

%% Comparison between all non-eHMI FAV scenarios (i.e motion-changing scenarios)
% This is comparing sc1, sc4 - sc7 (a total of 5 scenarios, in that order)
sc_idxs = [1,4,5,6,7];
mean_Likert = [];
mean_Mean_DTC = [];
mean_Cross_time = [];
mean_Gaze_ratio = [];
SEM_Likert = [];
SEM_Mean_DTC = [];
SEM_Cross_time = [];
SEM_Gaze_ratio = [];
for i = 1:length(sc_idxs)
    mean_Likert = [mean_Likert, mean(LIKERT_TRUST_RATINGS(:,sc_idxs(i)))];
    mean_Mean_DTC = [mean_Mean_DTC, mean(MEAN_DTU(:,sc_idxs(i)))];
    mean_Cross_time = [mean_Cross_time, mean(TIME_TO_CROSS(:,sc_idxs(i)))];
    mean_Gaze_ratio = [mean_Gaze_ratio, mean(GAZE_RATIO(:,sc_idxs(i)))];
    SEM_Likert = [SEM_Likert, SEM(LIKERT_TRUST_RATINGS(:,sc_idxs(i)))];
    SEM_Mean_DTC = [SEM_Mean_DTC,SEM(MEAN_DTU(:,sc_idxs(i)))];
    SEM_Cross_time = [SEM_Cross_time,SEM(TIME_TO_CROSS(:,sc_idxs(i)))];
    SEM_Gaze_ratio = [SEM_Gaze_ratio,SEM(GAZE_RATIO(:,sc_idxs(i)))];
end

% T-test
for i = 1:length(sc_idxs)
    for j = 1:length(sc_idxs)
        [H_LIKERT(i,j), P_LIKERT(i,j)]  = ttest2(LIKERT_TRUST_RATINGS(:,sc_idxs(i)), LIKERT_TRUST_RATINGS(:,sc_idxs(j)));
        [H_DTC(i,j), P_DTC(i,j)]  = ttest2(MEAN_DTU(:,sc_idxs(i)), MEAN_DTU(:,sc_idxs(j)));
        [H_TTC(i,j), P_TTC(i,j)]  = ttest2(TIME_TO_CROSS(:,sc_idxs(i)), TIME_TO_CROSS(:,sc_idxs(j)));
        [H_GAZE(i,j), P_GAZE(i,j)]  = ttest2(GAZE_RATIO(:,sc_idxs(i)), GAZE_RATIO(:,sc_idxs(j)));
    end
end
% Signifcance only for Likert and Gaze measures
fprintf("-------------------------------\n")
fprintf("Statistical Significance (Likert): sc1 vs sc4, p=%f\n",P_LIKERT(1,2))
fprintf("                                   sc1 vs sc6, p=%f (*)\n",P_LIKERT(1,4))
fprintf("                                   sc1 vs sc7, p=%f\n",P_LIKERT(1,5))
fprintf("                                   sc4 vs sc5, p=%f\n",P_LIKERT(2,3))
fprintf("                                   sc5 vs sc6, p=%f\n",P_LIKERT(3,4))
fprintf("                                   sc5 vs sc7, p=%f\n",P_LIKERT(3,5))

fprintf("Statistical Significance (Gaze):   sc1 vs sc6, p=%f (*)\n",P_GAZE(1,4))
fprintf("                                   sc4 vs sc6, p=%f\n",P_GAZE(2,4))
fprintf("                                   sc6 vs sc7, p=%f\n",P_GAZE(4,5))

% Plotting
y_labels = {"Likert Subjective-Trust Score (1-7)", "Mean DTC (uu)", "Cross Time (s)", "Gaze Ratio"};
x = 1:5;
x_labels = {"sc1", "sc4", "sc5", "sc6", "sc7"};
MEANS_combined = [mean_Likert;mean_Mean_DTC;mean_Cross_time;mean_Gaze_ratio];
SEMS_combined = [SEM_Likert;SEM_Mean_DTC;SEM_Cross_time;SEM_Gaze_ratio];
figure('Name','No eHMI FAV different Trajectories') 
set(gcf, 'Color', 'w'); set(gcf, 'Position',  [100, 100, 1100, 700]);
tiledlayout(2, 2)
for i = 1:4
    ax = nexttile();
    bar(ax, x,MEANS_combined(i,:)); hold on;
    e= errorbar(x,MEANS_combined(i,:), SEMS_combined(i,:)); grid on;e.Color = 'black'; e.LineWidth = 1;
    ylabel(y_labels{i}, 'FontSize',14);
    if i == 2
        ylim([220,260])
    end
    set(ax,'xticklabel',x_labels, 'fontsize',14); %xtickangle(gca,30);
end

% box and whisker plot
figure('Name','No eHMI FAV different Trajectories- Box plot') 
set(gcf, 'Color', 'w'); set(gcf, 'Position',  [100, 100, 1100, 700]);
combined_matrices(:,:,1) = LIKERT_TRUST_RATINGS;
combined_matrices(:,:,2) = MEAN_DTU;
combined_matrices(:,:,3) = TIME_TO_CROSS;
combined_matrices(:,:,4) = GAZE_RATIO;
tiledlayout(2, 2)
for i = 1:4
    ax = nexttile();
    boxplot(ax,combined_matrices(:,sc_idxs,i))
    ylabel(y_labels{i}, 'FontSize',14); grid on;
    set(ax,'xticklabel',x_labels, 'fontsize',14); %xtickangle(gca,30);
end

%% Connection to model parameters 
load("dict.mat");
traj_idxs = [8,11,14,29,55]; % The trajectories that were modelled in order [sc1,sc4,sc5,sc6,sc7]

% Obtain markers from trajectory dictionary (each marker in order of sc)
Vi = []; %initial velocity of the vehicle (km/hr)
Tbrake = []; %Time when brakes were first applied (s) - Note each trajectory starts from the same DTC (50m)
Dstop = []; %Distance to the crossing when vehicle stops (not relevant for crawl traj in sc7)
ACCmax = []; %Maximum deceleration (m/s^2)
FinalJerk = [];%(m/s^3) (not relevant to crawl traj in sc7 as final speed is constant)
MaxJerk = []; %("")
for i = 1:length(traj_idxs)
    Vi = [Vi,dict(traj_idxs(i)).v_brake];
    Tbrake = [Tbrake, dict(traj_idxs(i)).t_brake];
    Dstop = [Dstop, dict(traj_idxs(i)).d_stop];
    ACCmax = [ACCmax, dict(traj_idxs(i)).a_m];
    FinalJerk = [FinalJerk, dict(traj_idxs(i)).final_jerk];
    MaxJerk = [MaxJerk, dict(traj_idxs(i)).jerk_max];
end

% Print out scenario main vars
sc_labels = {"sc1", "sc4", "sc5", "sc6", "sc7"};
fprintf("\n--------------------------------------------------------------------------------\n")
fprintf("|  Sc  | v_i (km/hr) |  t_b (s) |  d_stop (m)  | v_f (km/hr) | a_max (m/s^2) |  \n")
fprintf("--------------------------------------------------------------------------------\n")
formatSpec = "|  %s  |      %i     |     %i    |       %i      |      %i      |      %.2f     | \n";
for i = 1:length(sc_labels)
    
    if i == 5
        vf = 5;
    else
        vf = 0;
    end

    fprintf(formatSpec, sc_labels{i}, Vi(i), Tbrake(i), Dstop(i), vf, ACCmax(i));
    fprintf("--------------------------------------------------------------------------------\n")
end
fprintf("--------------------------------------------------------------------------------\n")
fprintf("-------------------------Single Factor Changes Summary--------------------------\n")
fprintf("sc1 -> sc4: (*) distance to crossing at stop reduces by 2m\n")
fprintf("sc1 -> sc5: longer wait till brakes applied - slightly bigger decel\n")
fprintf("sc4 -> sc6: (*) higher V_i - bigger decel as a consequence\n")
fprintf("sc7 -> all: (***)crawling at 5km/hr\n")
fprintf("--------------------------------------------------------------------------------\n")


% Order x variables so that line plots are descriptive
[~, order_Vi]= sort(Vi);
[~,order_Tbrake] = sort(Tbrake);
[~,order_Dstop]= sort(Dstop);
[~,order_ACCmax]= sort(ACCmax);
[~,order_FinalJerk]= sort(FinalJerk);
[~,order_MaxJerk]= sort(MaxJerk);


% Get normalised partipant responses
sc_idxs = [1,4,5,6,7];
mean_nLikert = [];
mean_nMean_DTC = [];
mean_nCross_time = [];
mean_nGaze_ratio = [];
for i = 1:length(sc_idxs)
    mean_nLikert = [mean_nLikert, mean(norm_LIKERT_TRUST_RATINGS(:,sc_idxs(i)))];
    mean_nMean_DTC = [mean_nMean_DTC, mean(norm_MEAN_DTU(:,sc_idxs(i)))];
    mean_nCross_time = [mean_nCross_time, mean(norm_TIME_TO_CROSS(:,sc_idxs(i)))];
    mean_nGaze_ratio = [mean_nGaze_ratio, mean(norm_GAZE_RATIO(:,sc_idxs(i)))];
end

% Plotting relationships (normalised particpant markers so they can all be
% observed on the same axis)
figure('Name', 'Participant Responses and Trajectory Markers');set(gcf, 'Position',  [100, 100, 1100, 700]);
set(gcf,'Color','w')
orders = [order_Vi; order_Tbrake;order_Dstop;order_ACCmax;order_FinalJerk;order_MaxJerk];
x_vars = [Vi; Tbrake;Dstop;ACCmax;FinalJerk;MaxJerk];
x_labels = {"Initial Velocity (km/hr)", "Time when brakes are first applied (s)", ...
    "Final Distance to Crossing (m)","Maximum Acceleration (m/s^2)", "Final Jerk (m/s^3)", ...
    "Maximum Jerk (m/s^3)"}; 

for i = 1:length(x_vars)
    subplot(3,2,i)
    x = x_vars(i,:);

    if i == 3  %distance at stopping is irrelevant for sc7
        x(end) = [];
        [~, ord] = sort(x);
        x = x(ord);
        temp_mean_nLikert = mean_nLikert(1:end-1); %remove last data point
        temp_mean_nMean_DTC = mean_nMean_DTC(1:end-1);
        temp_mean_nCross_time =mean_nCross_time(1:end-1);
        temp_mean_nGaze_ratio = mean_nGaze_ratio(1:end-1);
        temp_mean_nLikert = temp_mean_nLikert(ord);
        temp_mean_nMean_DTC = temp_mean_nMean_DTC(ord); %order with new order
        temp_mean_nCross_time =temp_mean_nCross_time(ord);
        temp_mean_nGaze_ratio = temp_mean_nGaze_ratio(ord);
    else
        x = x(orders(i,:));
        temp_mean_nLikert = mean_nLikert(orders(i,:));
        temp_mean_nMean_DTC = mean_nMean_DTC(orders(i,:));
        temp_mean_nCross_time =mean_nCross_time(orders(i,:));
        temp_mean_nGaze_ratio = mean_nGaze_ratio(orders(i,:));
    end


    plot(x, temp_mean_nLikert, '-ok'); grid on; hold on;
    plot(x, temp_mean_nMean_DTC, '-go');
    plot(x, temp_mean_nCross_time, '-yo');
    plot(x, temp_mean_nGaze_ratio, '-bo');
    xlabel(x_labels{i}, 'FontSize',14)

    if i == 1 || i == 3 || i == 5
        ylabel("Mean Participant Marker (normalised)",'FontSize',14)
    end
    
    if i == 1
        legend("Likert","DTC", "TTC", "Gaze Ratio")
    end

end

%% Repeat above but remove sc7 
figure('Name', 'Participant Responses and Trajectory Markers - sc7 Removed');set(gcf, 'Position',  [100, 100, 1100, 700]);
set(gcf,'Color','w')
orders = [order_Vi; order_Tbrake;order_Dstop;order_ACCmax;order_FinalJerk;order_MaxJerk];
x_vars = [Vi; Tbrake;Dstop;ACCmax;FinalJerk;MaxJerk];
x_labels = {"Initial Velocity (km/hr)", "Time when brakes are first applied (s)", ...
    "Final Distance to Crossing (m)","Maximum Acceleration (m/s^2)", "Final Jerk (m/s^3)", ...
    "Maximum Jerk (m/s^3)"}; 

for i = 1:length(x_vars)
    subplot(3,2,i)
    x = x_vars(i,:);

    x(end) = [];
    [~, ord] = sort(x);
    x = x(ord);
    temp_mean_nLikert = mean_nLikert(1:end-1); %remove last data point
    temp_mean_nMean_DTC = mean_nMean_DTC(1:end-1);
    temp_mean_nCross_time =mean_nCross_time(1:end-1);
    temp_mean_nGaze_ratio = mean_nGaze_ratio(1:end-1);
    temp_mean_nLikert = temp_mean_nLikert(ord);
    temp_mean_nMean_DTC = temp_mean_nMean_DTC(ord); %order with new order
    temp_mean_nCross_time =temp_mean_nCross_time(ord);
    temp_mean_nGaze_ratio = temp_mean_nGaze_ratio(ord);


    plot(x, temp_mean_nLikert, '-ok'); grid on; hold on;
    plot(x, temp_mean_nMean_DTC, '-go');
    plot(x, temp_mean_nCross_time, '-mo');
    plot(x, temp_mean_nGaze_ratio, '-bo');
    xlabel(x_labels{i}, 'FontSize',14)

    if i == 1 || i == 3 || i == 5
        ylabel("Mean Participant Marker (normalised)",'FontSize',14)
    end
    
    if i == 1
        legend("Likert","DTC", "TTC", "Gaze Ratio")
    end

end


