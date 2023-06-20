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
% s = 1/6 - 2/((m+2)*(m+3))+ 1/((2*m+2)*(2*m+3));
% q = m^2/( (2*m+2)*(m+2) );

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
v_i_range = [10,15];
t1 = [0,1,2];
vf = vf; % take prev value, can change if needed

legend_labels = {};
plot_line_colors = {"r", "g", "b", "c", "k", "m"};
color_idx = 1;
DISTANCE = {};
VT = {};
for i = 1:length(v_i_range)
    
    for ii = 1:length(t1)
        % First calculate breaking time using v_i, t_1, constant v_f, d_max and
        % m
        td = double(t_d(v_i_range(i), vf, t1(ii), d_crossing, m_const));
    
        time_span = linspace(0, td, 30); % time variable
    
        % Calculate v_t based off t_d and v_i
        vt = double(v_t(v_i_range(i), vf, td, m_const,time_span));
        AOC = trapz(time_span, vt)/3.6
        p  = plot([0, t1(ii)], [v_i_range(i),v_i_range(i)]); % initial constant velocity
        p.Color = plot_line_colors{color_idx};

        time_span = time_span + t1(ii); % shift time range for intitial velocity curve
        p = plot(time_span, vt);
        p.Color = plot_line_colors{color_idx};
        color_idx = color_idx +1;
        legend_labels{end + 1} = strcat('v_i = ', num2str(v_i_range(i)), ...
            ', t_1 = ', num2str(t1(ii)), '=> t_d =  ', num2str(td));
    
        % distance-velocity info 
        time_span = linspace(0, td, 30); % reset timespan to start at 0
        dt = gradient(time_span); 
        distance = cumsum(vt.*dt).*10/36;
        distance_offset = v_i_range(i)*t1(ii);
        DISTANCE{end+1} = [0, distance_offset, distance+distance_offset];
        VT{end+1} = [v_i_range(i),v_i_range(i),vt];
    end

end
xlabel('Time (s)')
ylabel('Velocity (km/hr)')
title(strcat('v_f = ', num2str(vf), 'km/hr, d_{crossing} = ', ...
    num2str(d_crossing), 'm, m_{const} = ', num2str(m_const)), "FontSize", 16)
subtitle('Velocity with Time', 'FontSize', 14)
% legend(legend_labels)

subplot(1,2,2); hold on; grid on;
color_idx = 1; 
for i = 1:length(DISTANCE)
    d = DISTANCE{i};
    vt = VT{i};
    p = plot(d(1:2), vt(1:2)); % intital curve
    p.Color = plot_line_colors{color_idx};
    p = plot(d(3:end), vt(3:end));
    p.Color = plot_line_colors{color_idx};
    color_idx = color_idx +1;
end
xlabel('Distance (m)')
ylabel('Velocity (km/hr)')
title(strcat('v_f = ', num2str(vf), 'km/hr, d_{crossing} = ', ...
    num2str(d_crossing), 'm, m_{const} = ', num2str(m_const)), "FontSize", 16)
subtitle('Velocity with Distance', 'FontSize', 14)
% legend(legend_labels)



% ------------- 
% Decreasing t_d (for constant v_f, d_max) increases the sharpnest (or max
% acceleration) of the v_t curve.
% ---------------

%% Creating Dictionary of Literature markers from model

