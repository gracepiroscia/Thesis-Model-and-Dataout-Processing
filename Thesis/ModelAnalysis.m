% Author: Grace Piroscia
% 
% Parameter sweep of the deceleration model to determine the general
% behaviour and thus predicted outcomes for the piecewise model
%
clc;
clear;

%% Symbolic function(s)
% v_t(): output in km/hr
syms v_t(v_i, v_f, t_d, m, t)                    % r 
v_t(v_i, v_f, t_d, m, t) = v_i + 3.6*( ((1+2*m)^(2+1/m))/(4 * m^2) ) ...
    * (v_f - v_i)/(3.6 * t_d * ((1+2*m)^(2+1/m))/4 * 1/((2*m+2)*(m+2)) )  * ...% this term is am - max acceleration in m/s^2
    (t^2/t_d) * (0.5 - 2*(t/t_d)^m/(m+2)  + (t/t_d)^(2*m)/(2*m+2));

% Deceleration time - t_d - as a function of v_i, t_1 (t_1 = times
% when breaks are first applied), v_f, d_max (distance to crossing), and m
% t_d():   output in seconds (constants, s,q that depend on m)
syms t_d(v_i, v_f, t_1, d_max, m)
t_d(v_i, v_f, t_1, d_max, m) =  ( d_max - v_i*(5/18)*t_1 )/( (5/18)*v_i + ...
    ( ( 1/6 - 2/((m+2)*(m+3))+ 1/((2*m+2)*(2*m+3)) )/( m^2/((2*m+2)*(m+2)) ) )*(5/18)*(v_f - v_i) ); 
                                    % s/q
% Constants definition: 
% ---------------------
% s = 1/6 - 2/((m+2)*(m+3))+ 1/((2*m+2)*(2*m+3));
% q = m^2/( (2*m+2)*(m+2) );

% v_t for piecewise model
syms v_t_pw(v_i, v_f, t_decel, t_1, m, t)                    
v_t_pw(v_i, v_f, t_decel, t_1, m, t) = v_i + 3.6*( ((1+2*m)^(2+1/m))/(4 * m^2) ) ...
    * (v_f - v_i)/(3.6 * t_decel * ((1+2*m)^(2+1/m))/4 * 1/((2*m+2)*(m+2)) )  * ...
    ((t-t_1)^2/t_decel) * (0.5 - 2*((t-t_1)/t_decel)^m/(m+2)  + ((t-t_1)/t_decel)^(2*m)/(2*m+2));
%% Validate with existing data
% % (from DecelModelFitting.m plots)
% m_constant = 5.5; 
% vf = 0.45; %(km/hr)
% vi = 19.5; %""
% td = 6.3; %(s)
% tVals = 0:0.01:td;
% 
% vVals = double(v_t(vi,vf,td, m_constant, tVals));
% plot(tVals, vVals)
% 
% d_travelled = trapz(tVals, vVals); %units m/3.6


%% Param sweep:
%% Constant final velocity (v_f) and distance to crossing (d_max)
d_crossing = 35; % Distance to Crossing (m)
vf = 5; % Final Velocity (km/hr)
m_range = linspace(0.05,19.05,5); % Range of m fitting parameters. Range determined from distribution of
                        % fit data found in DecelModelFitting.m, along with
                        % minimum found from Akcelik and Biggs (1987,
                        % Vol21, pp 36-54).
%% Effect on t_d:
% Adjusting v_i only:
v_i_range = linspace(2,60, 30); %km/hr
t1 = 2; % Time when brakes are first applied (s)

figure(1); set(gcf,'Color','w'); subplot (1,3,1); hold on; grid on;
legend_labels = {};
for i = 1:length(m_range)
    td = double(t_d(v_i_range, vf, t1, d_crossing, m_range(i)));
    plot(v_i_range, td);
    legend_labels{end+1} = strcat('m = ', num2str(m_range(i)));
end 
xlabel('Initial Velocity (km/hr)')
ylabel('Deceleration Time (s)')
title(strcat('v_f = ', num2str(vf), 'km/hr, t_{brake} = ', num2str(t1), 's , d_{crossing} = ', ...
    num2str(d_crossing), 'm '), "FontSize", 16)
legend(legend_labels)

% Adjusting t_1 only:
t_1_range = linspace(0, 5, 30);
vi = 40; % constant initial velocity 
subplot(1,3,2); hold on; grid on;
for i = 1:length(m_range)
    td = double(t_d(vi, vf, t_1_range, d_crossing, m_range(i)));
    plot(t_1_range, td);
end 
xlabel('Time when brakes are first applied (s)')
ylabel('Deceleration Time (s)')
title(strcat('v_f = ', num2str(vf), 'km/hr, v_i = ', num2str(vi), 'km/hr , d_{crossing} = ', ...
    num2str(d_crossing), 'm '), "FontSize", 16)
legend(legend_labels)

% Adjusting both v_i and t_1:
subplot(1,3,3)
for i = 1:length(m_range)
    td = double(t_d(v_i_range, vf, t_1_range, d_crossing, m_range(i)));
    plot3(v_i_range, t_1_range, td, '.-');
    hold on;
    drawnow
end 
grid on;
xlabel('Initial Velocity (km/hr)')
ylabel('Time when brakes are first applied (s)')
zlabel('Deceleration Time (s)')
title(strcat('v_f = ', num2str(vf), 'km/hr, d_{crossing} = ', ...
    num2str(d_crossing), 'm '), "FontSize", 16)
legend(legend_labels)

% ------------------------------------------------------------------------
% Increasing either v_i, t_1 or both, will result in a decrease in t_d for 
% constant final velocity, over a constant distance travelled
% ------------------------------------------------------------------------
%% Effect on v_t:
% Adjusting t_d and v_i only:
figure(2); set(gcf,'Color','w'); subplot(1,2,1); grid on; hold on;
m_const = 5; % convergent behaviour for m > 4 as seen in the above plots
v_i_range = linspace(10,40,6);
t1 = 0; % brakes applied straight away
vf = vf; % take prev value, can change if needed

legend_labels = {};
DISTANCE = {};
VT = {};
for i = 1:length(v_i_range)
    % First calculate breaking time using v_i, constant v_f, t_1, d_max and
    % m
    td = double(t_d(v_i_range(i), vf, t1, d_crossing, m_const));

    time_span = linspace(0, td, 30); % time variable

    % Calculate v_t based off t_d and v_i
    vt = double(v_t(v_i_range(i), vf, td, m_const,time_span));
    plot(time_span, vt)
    legend_labels{end + 1} = strcat('v_i = ', num2str(v_i_range(i)), '=> t_d =  ', num2str(td));

    % distance-velocity info 
    dt = gradient(time_span); 
    distance = cumsum(vt.*dt).*10/36;
    DISTANCE{end+1} = distance;
    VT{end+1} = vt;

end
xlabel('Time (s)')
ylabel('Velocity (km/hr)')
title(strcat('v_f = ', num2str(vf), 'km/hr, d_{crossing} = ', ...
    num2str(d_crossing), 'm, m_{const} = ', num2str(m_const)), "FontSize", 16)
subtitle('Velocity with Time', 'FontSize', 14)
legend(legend_labels)

subplot(1,2,2); hold on; grid on;
for i = 1:length(DISTANCE)
    plot(DISTANCE{i}, VT{i})
end
xlabel('Distance (m)')
ylabel('Velocity (km/hr)')
title(strcat('v_f = ', num2str(vf), 'km/hr, d_{crossing} = ', ...
    num2str(d_crossing), 'm, m_{const} = ', num2str(m_const)), "FontSize", 16)
subtitle('Velocity with Distance', 'FontSize', 14)
legend(legend_labels)



% ------------- 
% Decreasing t_d (for constant v_f, d_max) increases the sharpnest (or max
% acceleration) of the v_t curve.
% ---------------

%% Effect on v_t - introducing piecewise plots
figure(3); set(gcf,'Color','w'); subplot(1,2,1); grid on; hold on;
m_const = 5; % convergent behaviour for m > 4 as seen in the above plots
v_i_range = [10,15]; %km/hr
t1 = [0,2,4];  %(s) 
vf = vf; % take prev value, can change if needed

legend_labels = {};
plot_line_colors = {"r", "g", "b", "c", "k", "m"};
color_idx = 1;
DISTANCE = {};
VT = {};
for i = 1:length(v_i_range)
    
    for ii = 1:length(t1)
        % Constant velocity piece:
        p  = plot([0, t1(ii)], [v_i_range(i),v_i_range(i)]); % initial constant velocity
        p.Color = plot_line_colors{color_idx};        


        % For deceleration piece, first calculate breaking time using v_i, t_1, constant v_f, d_max and
        % m
        td = double(t_d(v_i_range(i), vf, t1(ii), d_crossing, m_const));
        time_span = linspace(t1(ii), td + t1(ii), 30); % time variable
        vt = double(v_t_pw(v_i_range(i), vf, td, t1(ii), m_const,time_span)); % Calculate v_t based off t_d and v_i
        p = plot(time_span, vt);
        p.Color = plot_line_colors{color_idx};
        color_idx = color_idx +1;
        legend_labels{end + 1} = strcat('v_i = ', num2str(v_i_range(i)), ...
            ', t_1 = ', num2str(t1(ii)), '=> t_d =  ', num2str(td));
    
        % distance-velocity info 
        dt = gradient(time_span); 
        distance = cumsum(vt.*dt).*10/36;
        distance_offset = v_i_range(i)*5/18*t1(ii); %(m) for distance travelled initially
        DISTANCE{end+1} = [0, distance_offset, distance+distance_offset];
        VT{end+1} = [v_i_range(i),v_i_range(i),vt];
    end

end
xlabel('Time (s)')
ylabel('Velocity (km/hr)')
title(strcat('v_f = ', num2str(vf), 'km/hr, d_{crossing} = ', ...
    num2str(d_crossing), 'm, m_{const} = ', num2str(m_const)), "FontSize", 16)
subtitle('Velocity with Time', 'FontSize', 14)

subplot(1,2,2); hold on; grid on;
color_idx = 1; 
for i = 1:length(DISTANCE)
    p = plot(DISTANCE{i}, VT{i});
    p.Color = plot_line_colors{color_idx};
    color_idx = color_idx +1;
end
xlabel('Distance (m)')
ylabel('Velocity (km/hr)')
title(strcat('v_f = ', num2str(vf), 'km/hr, d_{crossing} = ', ...
    num2str(d_crossing), 'm, m_{const} = ', num2str(m_const)), "FontSize", 16)
subtitle('Velocity with Distance', 'FontSize', 14)
legend(legend_labels)


%% Creating Dictionary of Literature markers from model
% Velocity when brakes are first applied - v_brake (km/hr)
% Time when brakes are first applied - t_brake (s)
% Calibration parameter m 
% Deceleration time (t_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Max acceleration - a_max (m/s^2)
% Maximum jerk - jerk_max (m/s^3)
% Stopping distance (final dist between pedestrian and vehicle) - d_stop (m)
% Final velocity - v_final (m/s)
% Distance-to-crossing-brake (distance to crossing when brakes are first applied) - DTC_brake
% Time-to-zebra-brake (time to pedestrian crossing when brakes are first applied) - TTZ_brake (s)

d_range = 50;% All simulations run over the same range from crossing(m)
d_traj = [47,49]; % distance travelled by vehicle (should always be < d_range as greater would indicate pedestrian collision) (m)
v_i_range = [30,40,50]; %km/hr
t1 = [0,2,3];  %Time when brakes are first applied (s) 
vf = 0; % All stopping trajectories 
m_range = linspace(0.05,19.05,3);

i_dict = 1;
for i = 1:length(v_i_range)
    
    for ii = 1:length(t1)
        for iii = 1:length(d_traj)
            for iiii = 1:length(m_range)
       
                % For deceleration piece, first calculate breaking time using v_i, t_1, m, constant v_f and d_max
                td = double(t_d(v_i_range(i), vf, t1(ii), d_traj(iii), m_range(iiii)));
                time_span = linspace(t1(ii), td + t1(ii), 40); % time variable
    
                dict(i_dict).v_brake = v_i_range(i); % Store
                dict(i_dict).t_brake = t1(ii);
                dict(i_dict).tot_dist_travelled = d_traj(iii);
                dict(i_dict).m = m_range(iiii);
                dict(i_dict).t_d = td;
                dict(i_dict).a_m = (vf - v_i_range(i))/(3.6 * td * ...% using relation instead of measuring directly
                    ((1+2*m_range(iiii))^(2+1/m_range(iiii)))/4 * 1/((2*m_range(iiii)+2)*(m_range(iiii)+2)) );
                dict(i_dict).a_m = abs(dict(i_dict).a_m);
    
                % Now calculate v_t based off t_d and v_i 
                vt = double(v_t_pw(v_i_range(i), vf, td, t1(ii), m_range(iiii),time_span)); 
    
                % Calculate acceleration series using vt
                dt = gradient(time_span);
                at = gradient(vt)./dt;
                jerk = gradient(at)./dt(1);
                dict(i_dict).jerk_max = max(abs(jerk));
                dict(i_dict).final_jerk = abs(jerk(end));
                dict(i_dict).d_stop = d_range - d_traj(iii);
                dict(i_dict).v_final = vf;
                dict(i_dict).DTC_brake = d_range - v_i_range(i)*5/18*t1(ii); %range minus distance travelled during constant traj
                dict(i_dict).TTZ_brake = dict(i_dict).DTC_brake/v_i_range(i);

                % Store trajectory
                dict(i_dict).vt_traj = [v_i_range(i),v_i_range(i),vt];
                dict(i_dict).time = [0,t1(ii), time_span];
                distance = cumsum(vt.*dt).*10/36;
                distance_offset = v_i_range(i)*5/18*t1(ii); %(m) for distance travelled initially
                dict(i_dict).vd_traj = [0, distance_offset, distance+distance_offset];

                % In a reasonable acceleration range:
                if (dict(i_dict).a_m >= 2) && (dict(i_dict).a_m <= 5)
                    dict(i_dict).is_reasonable_traj = 1;
                else
                    dict(i_dict).is_reasonable_traj = 0;
                end
                
                i_dict = i_dict +1;
    
            end 
        end
    end

end
%1 crawl scenatio with v_i = 30, t_b = 2, d_to_crossing_when_crawl_begins
%=4m, v_f = 5km/hr
vi_crawl = 30;
t_b_crawl = 2;
vf_crawl = 5;
d_traj_crawl = d_range - 5;
m_crawl = 9.55;


% For deceleration piece, first calculate breaking time using v_i, t_1, m, constant v_f and d_max
td = double(t_d(vi_crawl, vf_crawl, t_b_crawl, d_traj_crawl, m_crawl));
time_span = linspace(t_b_crawl, td + t_b_crawl, 40); % time variable

dict(i_dict).v_brake = vi_crawl; % Store
dict(i_dict).t_brake = t_b_crawl;
dict(i_dict).tot_dist_travelled = d_traj_crawl;
dict(i_dict).m = m_crawl;
dict(i_dict).t_d = td;
dict(i_dict).a_m = (vf_crawl - vi_crawl)/(3.6 * td * ...% using relation instead of measuring directly
    ((1+2*m_crawl)^(2+1/m_crawl))/4 * 1/((2*m_crawl+2)*(m_crawl+2)) );
dict(i_dict).a_m = abs(dict(i_dict).a_m);

% Now calculate v_t based off t_d and v_i 
vt = double(v_t_pw(vi_crawl, vf_crawl, td, t_b_crawl, m_crawl,time_span)); 

% Calculate acceleration series using vt
dt = gradient(time_span);
at = gradient(vt)./dt;
jerk = gradient(at)./dt(1);
dict(i_dict).jerk_max = max(abs(jerk));
dict(i_dict).final_jerk = abs(jerk(end));
dict(i_dict).d_stop = d_range - d_traj_crawl;
dict(i_dict).v_final = vf_crawl;
dict(i_dict).DTC_brake = d_range - vi_crawl*5/18*t_b_crawl; %range minus distance travelled during constant traj
dict(i_dict).TTZ_brake = dict(i_dict).DTC_brake/vi_crawl;

% Store trajectory
dict(i_dict).vt_traj = [vi_crawl,vi_crawl,vt];
dict(i_dict).time = [0,t_b_crawl, time_span];
distance = cumsum(vt.*dt).*10/36;
distance_offset = vi_crawl*5/18*t_b_crawl; %(m) for distance travelled initially
dict(i_dict).vd_traj = [0, distance_offset, distance+distance_offset];

% In a reasonable acceleration range:
if (dict(i_dict).a_m >= 2) && (dict(i_dict).a_m <= 5)
    dict(i_dict).is_reasonable_traj = 1;
else
    dict(i_dict).is_reasonable_traj = 0;
end


save("dict.mat", "dict");


