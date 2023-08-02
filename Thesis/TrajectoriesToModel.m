% Author: Grace Piroscia
% Plot the trajectories and distance-frames information for modelling in
% maya.
%
% Following helpful definitions:
% - v_i: Initial velocity of the vehicle (km/hr)
% - t_b: Time in the trajectory when brakes a first applied (s)
% - d_stop: The distance between the pedestrian crossing and the car when
%   the vehicle has reached a stop (m)
% - v_f: The final velocity of the vehicle (km/hr)
% - a_m: the maximum deceleration achieved in the vehicles trajecotry
%   (m/s^2)
% - m: model parameter
%
clc;
clear;

%% Load data 
load("dict.mat")

%% Idx's for ideal trajectories
idxs = [8,11,14,29,55];

%% Lets print out the data to be clear the difference in trajectories
fprintf("-------------------------------------------------------------------------------------\n")
fprintf("|  #  | v_i (km/hr) |  t_b (s) |  d_stop (m)  | v_f (km/hr) | a_max (m/s^2) |    m    \n")
fprintf("-------------------------------------------------------------------------------------\n")
formatSpec = "|  %i  |      %i     |     %i    |       %i      |      %i      |      %.2f     |  %.2f  \n";
for i = 1:length(idxs)
    
    vi = dict(idxs(i)).v_brake;
    tb = dict(idxs(i)).t_brake;
    dstop = dict(idxs(i)).d_stop;
    vf = dict(idxs(i)).v_final;
    am = dict(idxs(i)).a_m;
    m_ = dict(idxs(i)).m;

    fprintf(formatSpec, i, vi, tb, dstop, vf, am, m_);
    fprintf("-------------------------------------------------------------------------------------\n")
end

%% Velocity-distance plots
% legend lables format
llabels = {};
for i = 1:length(idxs)
    str_lable = "idx = " + num2str(idxs(i));
    llabels{i} = str_lable;
end
figure(1)
set(gcf, "Color", 'w')
hold on
for i = 1:length(idxs)
    %plot velocity with distance:
    dist = dict(idxs(i)).vd_traj;
    vel  = dict(idxs(i)).vt_traj;

    plot(dist,vel)

end
xlabel("Distance (m)")
ylabel("Velocity (km/hr)")
legend(llabels)
grid on;

%% Distance-time plots
figure(2)
set(gcf, "Color", 'w')
hold on
for i = 1:length(idxs)
    %plot velocity with distance:
    dist = dict(idxs(i)).vd_traj;
    time_s  = dict(idxs(i)).time;

    plot(time_s,dist)

end
xlabel("time (s)")
ylabel("Distance (m)")
legend(llabels)
grid on;

%% Distance -frames plots for help with maya construction
figure(3)
set(gcf, "Color", 'w')
hold on
for i = 1:length(idxs)
    %plot velocity with distance:
    dist = dict(idxs(i)).vd_traj; 
    time_s  = dict(idxs(i)).time;
    frames = time_s.*24;

    plot(frames,dist)

end
xlabel("Frames (24fps)")
ylabel("Distance (m)")
legend(llabels)
grid on;