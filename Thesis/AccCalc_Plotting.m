clc;
clear;

load("Combined_Trajectories.mat");


%% Calculating Acceleration and distance
%Crawl Trajectories
for i = 1:length(combined.CrawlTrajectories)
    
    VX = combined.CrawlTrajectories(i).Vx;
    t = combined.CrawlTrajectories(i).time_s;
    
    % Acceleration
    Ax = gradient(VX)./gradient(t);
    combined.CrawlTrajectories(i).Ax_ms2 = Ax;
    
    % Distance travelled until crossing (taken as the distance at the final velocity value)
    dt = gradient(t);
    dS = VX.*dt;  %element-wise
    d = cumsum(dS);     %cumulative distance
    
    combined.CrawlTrajectories(i).distance_m = d;
    
    
end

% Stop Trajectories
for i = 1:length(combined.StopTrajectories)
    
    VX = combined.StopTrajectories(i).Vx;
    t = combined.StopTrajectories(i).time_s;
    
    % Acceleration
    Ax = gradient(VX)./gradient(t);
    combined.StopTrajectories(i).Ax_ms2 = Ax;
    
    % Distance travelled until crossing (taken as the distance at the final velocity value)
    dt = gradient(t);
    dS = VX.*dt;  %element-wise
    d = cumsum(dS);     %cumulative distance
    combined.StopTrajectories(i).distance_m = d;
    
    
end

% save struct
%save("Combined_Trajectories.mat", 'combined')


%% Plotting velocity vs time:
% Crawl
figure(1)
subplot(2,1,1) %single pedestrian
hold on
for i = 1:length(combined.CrawlTrajectories)
    
    VX = combined.CrawlTrajectories(i).Vx;
    time = combined.CrawlTrajectories(i).time_s;
    
    if combined.CrawlTrajectories(i).isSinglePedestrian == 1
        plot(time, VX,'r')
    else
        plot(time, VX, 'k')
    end
        
end
grid on;
set(gcf, 'color', 'w')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Crawl Trajectories')


% Stop trajectories
subplot(2,1,2)
hold on
for i = 1:length(combined.StopTrajectories)
    
    VX = combined.StopTrajectories(i).Vx;
    time = combined.StopTrajectories(i).time_s;
    
    if combined.StopTrajectories(i).isSinglePedestrian == 1
        plot(time, VX,'r')
    else
        plot(time, VX, 'k')
    end
        
end
grid on;
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Stop Trajectories')

%% Plotting: Velocity VS Distance
% Crawl
figure(2)
subplot(2,1,1) %single pedestrian
hold on
for i = 1:length(combined.CrawlTrajectories)
    
    VX = combined.CrawlTrajectories(i).Vx;
    Distance = combined.CrawlTrajectories(i).distance_m;
    
    if combined.CrawlTrajectories(i).isSinglePedestrian == 1
        plot(Distance, VX,'r')
    else
        plot(Distance, VX, 'k')
    end
        
end
grid on;
set(gcf, 'color', 'w')
xlabel('Distance (m)')
ylabel('Velocity (m/s)')
title('Crawl Trajectories')


% Stop trajectories
subplot(2,1,2)
hold on
for i = 1:length(combined.StopTrajectories)
    
    VX = combined.StopTrajectories(i).Vx;
    Distance = combined.StopTrajectories(i).distance_m;
    
    if combined.StopTrajectories(i).isSinglePedestrian == 1
        plot(Distance, VX,'r')
    else
        plot(Distance, VX, 'k')
    end
        
end
grid on;
xlabel('Distance (m)')
ylabel('Velocity (m/s)')
title('Stop Trajectories')

%% Plotting: Acceleration VS Distance
% Crawl
figure(3)
subplot(2,1,1) %single pedestrian
hold on
for i = 1:length(combined.CrawlTrajectories)
    
    AX = combined.CrawlTrajectories(i).Ax_ms2;
    Distance = combined.CrawlTrajectories(i).distance_m;
    
    if combined.CrawlTrajectories(i).isSinglePedestrian == 1
        plot(Distance, AX,'r')
    else
        plot(Distance, AX, 'k')
    end
        
end
grid on;
set(gcf, 'color', 'w')
xlabel('Distance (m)')
ylabel('Acceleration (m/s^2)')
title('Crawl Trajectories')


% Stop trajectories
subplot(2,1,2)
hold on
for i = 1:length(combined.StopTrajectories)
    
    AX = combined.StopTrajectories(i).Ax_ms2;
    Distance = combined.StopTrajectories(i).distance_m;
    
    if combined.StopTrajectories(i).isSinglePedestrian == 1
        plot(Distance, AX,'r')
    else
        plot(Distance, AX, 'k')
    end
        
end
grid on;
xlabel('Distance (m)')
ylabel('Acceleration (m/s^2)')
title('Stop Trajectories')
