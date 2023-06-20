% Author: Grace Piroscia
%
% Scipt to fit deceleration model to EV data
% understanding how parameters interact dynamically with this model
%
% Things to note so far:
% - Model below is calibrated for real-life vehicles at a signalised intersection
%
%TODO:
% 0. What to do with nudging trajectories? Do we just say these aren't
% accounted for in the model and a more complex trajectory could be used (i.e. just more of a 
% piece-wise analysis at each deceleration event) I think this would be best.

clc
clf;
clear

%% Load relevant data
load("Combined_Trajectories.mat");


%% Equation - velocity with time

%v_t = 'v_i + 3.6*( ((1+2*m)^(2+1/m))/(4 * m^2) ) * a_m  * (t^2/t_d) * (0.5 - 2*(t/t_d)^m/(m+2)  + (t/t_d)^(2*m)/(2*m+2))';

v_t = 'v_i + 3.6*( ((1+2*m)^(2+1/m))/(4 * m^2) ) * (v_f - v_i)/(3.6*t_d * ((1+2*m)^(2+1/m))/4 * 1/((2*m+2)*(m+2)) )  * (t^2/t_d) * (0.5 - 2*(t/t_d)^m/(m+2)  + (t/t_d)^(2*m)/(2*m+2))';


%% Fitting
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

% range for parameter m
% opts.Lower = [0,0];
% opts.StartPoint = [1,1];
% opts.upper = [inf,inf];
opts.Lower = 0;
opts.StartPoint = 1;
opts.upper = inf;


% range for parameters a_m, m, t_d (alphabetised)
% opts.Lower = [-inf,0,0];
% opts.StartPoint = [-1,1,2];
% opts.upper = [0,inf,inf];

fn = fittype(v_t, 'independent', 't', 'problem', {'v_i', 'v_f', 't_d'}, 'options', opts);
%fn = fittype(v_t, 'independent', 't', 'options', opts);

for i = 1:length(combined.StopTrajectories)
    %% Constants 
    % From EV data
    v_data = combined.StopTrajectories(i).Vx * 3.6; %(km/h)
    t_data = combined.StopTrajectories(i).time_s;
    acc_data = combined.StopTrajectories(i).Ax_ms2; %(m/s^2)
    d_data = combined.StopTrajectories(i).distance_m;% (m)
    
    [acc_data, v_data, t_data, d_data] = ProcessTrajectory(acc_data, v_data, t_data, d_data);
   
    t_i = t_data(1); %(s)
    t_f = t_data(end); %""
    v_i = v_data(1); % (km/hr)
    v_f = v_data(end); %""
    t_d = t_f - t_i;    %we can directly measure this
    a_m = -max(abs(acc_data));%-1.5; %max deceleration rate (m/s^2)
     

    % fit;
    [fitresult, gof] = fit(t_data,v_data,fn, 'problem', {v_i, v_f, t_d} );
    coefficientValues = coeffvalues(fitresult);
    %[fitresult, gof] = fit(t_data,v_data,fn);
    
    % for comparison to measured val
    m = coefficientValues;
    a_m_calculated = (v_f - v_i)/(3.6*t_d * ((1+2*m)^(2+1/m))/4 * 1/((2*m+2)*(m+2)) );
    
    time_step = t_i:0.05:t_f; %(s)
    dt = gradient(time_step); 
    fitted_velocity = feval(fitresult, time_step);
    fitted_distance = cumsum(fitted_velocity.*dt').*10/36;

    %% Plotting
    % figure(i)
    % set(gcf,'color', 'w');
    % plot(fitresult); grid on;
    % hold on
    % plot(t_data,v_data, 'm');
    % xlabel('Time (s)')
    % ylabel('Velocity (km/hr)')
    % title(["i = " + num2str(i)])
    % hold off
    % 
    % figure(i+1)
    % hold on
    % set(gcf,'color', 'w');
    % plot(fitted_distance,fitted_velocity, 'b'); grid on; hold on; 
    % plot(d_data, v_data, 'm');
    % xlabel('Distance (m)')
    % ylabel('Velocity (km/hr)')
    % title(["i = " + num2str(i)])
    % hold off

    % store data
    fit_data.StopTraject(i).vi_km_h = v_i;
    fit_data.StopTraject(i).am_m_s_s = a_m;
    fit_data.StopTraject(i).td_s = t_d;
    fit_data.StopTraject(i).m_fit = coefficientValues; % 95% confidence bounds
    fit_data.StopTraject(i).std_error = gof.rmse; %root mean square error (or std err)
    fit_data.StopTraject(i).r_sqr = gof.rsquare;


end 

% Repeat for crawl trajectories
for i = 1:length(combined.CrawlTrajectories)
    %% Constants 
    % From EV data
    v_data = combined.CrawlTrajectories(i).Vx * 3.6; %(km/h)
    t_data = combined.CrawlTrajectories(i).time_s;
    acc_data = combined.CrawlTrajectories(i).Ax_ms2; %(m/s^2)
    d_data = combined.CrawlTrajectories(i).distance_m;% (m)
    
    [acc_data, v_data, t_data, d_data] = ProcessTrajectory(acc_data, v_data, t_data, d_data);
   
    t_i = t_data(1); %(s)
    t_f = t_data(end); %""
    v_i = v_data(1); % (km/hr)
    v_f = v_data(end); %""
    a_m = -max(abs(acc_data));%-1.5; %max deceleration rate (m/s^2)
    t_d = t_f - t_i;

    % fit
    [fitresult, gof] = fit(t_data,v_data,fn, 'problem', {v_i, v_f, t_d} );
    coefficientValues = coeffvalues(fitresult);
    
    time_step = t_i:0.05:t_f; %(s)
    dt = gradient(time_step); 
    fitted_velocity = feval(fitresult, time_step);
    fitted_distance = cumsum(fitted_velocity.*dt').*10/36;

    %% Plotting
    % figure(i)
    % set(gcf,'color', 'w');
    % plot(fitresult); grid on;
    % hold on
    % plot(t_data,v_data, 'm');
    % xlabel('Time (s)')
    % ylabel('Velocity (km/hr)')
    % title(["i = " + num2str(i)])
    % hold off

    % store data
    fit_data.CrawlTraject(i).vi_km_h = v_i;
    fit_data.CrawlTraject(i).am_m_s_s = a_m;
    fit_data.CrawlTraject(i).td_s = t_d;
    fit_data.CrawlTraject(i).m_fit = coefficientValues; % 95% confidence bounds
    fit_data.CrawlTraject(i).std_error = gof.rmse; %root mean square error (or std err)
    fit_data.CrawlTraject(i).r_sqr = gof.rsquare;


end 
%% Plot Distribution of m-fit values and errors
StopTraject_m = vertcat(fit_data.StopTraject.m_fit);
CrawlTraject_m = vertcat(fit_data.CrawlTraject.m_fit);
CrawlTraject_m(11)= []; %rm outlier
figure(1)
set(gcf,'color', 'w');
subplot(1,2,1)
bar(1:length(StopTraject_m), StopTraject_m); grid on;% histogram(StopTraject_m); grid on; 
title('Stop Trajectories')
ylabel('m')
xlabel('idx')
subplot(1,2,2)
bar(1:length(CrawlTraject_m), CrawlTraject_m); grid on;%histogram(CrawlTraject_m); grid on;
title('Crawl Trajectories')
ylabel('m')
xlabel('idx')

StopTraject_rsqr = vertcat(fit_data.StopTraject.r_sqr);
CrawlTraject_rsqr = vertcat(fit_data.CrawlTraject.r_sqr);
CrawlTraject_rsqr(11)= []; %rm outlier
figure(2)
set(gcf,'color', 'w');
subplot(1,2,1)
bar(1:length(StopTraject_rsqr), StopTraject_rsqr);grid on;%histogram(StopTraject_rsqr); grid on;
title('Stop Trajectories')
ylabel('Rsqr')
xlabel('idx')
subplot(1,2,2)
bar(1:length(CrawlTraject_rsqr), CrawlTraject_rsqr);grid on;%histogram(CrawlTraject_rsqr); grid on;
title('Crawl Trajectories')
ylabel('Rsqr')
xlabel('idx')

%% Narrow down values
r_sqr_thresh = 0.9;
Narrowed_stopTraject_rsqr = StopTraject_rsqr(StopTraject_rsqr >=r_sqr_thresh);
Narrowed_stopTraject_m = StopTraject_m(StopTraject_rsqr >=r_sqr_thresh);

Narrowed_crawlTraject_rsqr = CrawlTraject_rsqr(CrawlTraject_rsqr >=r_sqr_thresh);
Narrowed_crawlTraject_m = CrawlTraject_m(CrawlTraject_rsqr >=r_sqr_thresh);
figure(3)
set(gcf,'color', 'w');
subplot(1,2,1)
bar(1:length(Narrowed_stopTraject_m), Narrowed_stopTraject_m); grid on;% histogram(StopTraject_m); grid on; 
title('Stop Trajectories')
subtitle(strcat('RSQR >= ',num2str(r_sqr_thresh)))
ylabel('m')
xlabel('idx')
subplot(1,2,2)
bar(1:length(Narrowed_crawlTraject_m), Narrowed_crawlTraject_m); grid on;%histogram(CrawlTraject_m); grid on;
title('Crawl Trajectories')
subtitle(strcat('RSQR >= ',num2str(r_sqr_thresh)))
ylabel('m')
xlabel('idx')



