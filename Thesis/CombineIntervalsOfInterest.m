% Author: Grace Piroscia
%
% Combine 'crawl' and 'stop' trajectories from multiple different bags into
% a single struct for later plotting and fitting. Note: this requires the
% time interval of interests to be known already, as done with
% TimeIntervalsOfInterest.m script.

clc;
clear;

FileNames = ["2018-03-08-14-30-07_Dataset_year_.mat",...
    "2018-03-16-10-33-39_Dataset_year.mat",...
    "2018-03-22-10-08-56_Dataset_year.mat",...
    "2018-03-27-09-38-01_Dataset_Year.mat",...
    "2018-04-05-13-28-38_Dataset_week.mat",...
    "2018-04-10-15-45-04_Dataset_year.mat",...
    "2018-04-19-15-05-55_Dataset_Year.mat"];

i_crawl = 1; %struct indices
i_stop = 1;

for File = 1:length(FileNames)
    
    FileName = FileNames(File);
    load(FileName);

    VX = vertcat(s.odomRear.LinearVel_X);
    time = vertcat(s.odomRear.Time_since_start_secs);
    time = round(time,2); %round to 2dp to for exact find of time of interest

    % time intervals for trajectories of interest only
    timeIntervalOfInterest_crawl = s.odomRear.Crossing_Crawl_event_timeSinceStart;
    timeIntervalOfInterest_stop =  s.odomRear.Crossing_Complete_Stop_event_timeSinceStart;
    
    isSinglePedestrian_stop = s.odomRear.IsSinglePedestrianStop;
    isSinglePedestrian_crawl = s.odomRear.IsSinglePedestrianStopCrawl;

    [~, ncols] = size(timeIntervalOfInterest_crawl);
    for i = 1:ncols

        start_time = timeIntervalOfInterest_crawl(1,i);
        end_time = timeIntervalOfInterest_crawl(2,i);

        start_idx = find(time == start_time);  
        end_idx = find(time == end_time);

        cut_vx = VX(start_idx:end_idx);
        cut_time = time(start_idx:end_idx);
        cut_time_relative = cut_time - cut_time(1); % time relative to start
        
        combined.CrawlTrajectories(i_crawl).Vx = cut_vx;
        combined.CrawlTrajectories(i_crawl).time_s = cut_time_relative;
        
        if isSinglePedestrian_crawl(i) == 1
            combined.CrawlTrajectories(i_crawl).isSinglePedestrian = 1;
        else
            combined.CrawlTrajectories(i_crawl).isSinglePedestrian = 0;
        end
        
        combined.CrawlTrajectories(i_crawl).source = FileName;
        
        i_crawl = i_crawl + 1;
        
    end

    [~, ncols] = size(timeIntervalOfInterest_stop);
    for i = 1:ncols

        start_time = timeIntervalOfInterest_stop(1,i);
        end_time = timeIntervalOfInterest_stop(2,i);

        start_idx = find(time == start_time);  
        end_idx = find(time == end_time);

        cut_vx = VX(start_idx:end_idx);
        cut_time = time(start_idx:end_idx);
        cut_time_relative = cut_time - cut_time(1); % time relative to start
        
        combined.StopTrajectories(i_stop).Vx = cut_vx;
        combined.StopTrajectories(i_stop).time_s = cut_time_relative;
        
        if isSinglePedestrian_stop(i) == 1
            combined.StopTrajectories(i_stop).isSinglePedestrian = 1;
        else
            combined.StopTrajectories(i_stop).isSinglePedestrian = 0;
        end
        
        combined.StopTrajectories(i_stop).source = FileName;
        
        i_stop = i_stop + 1;

    end

end 

%% Store
save("Combined_Trajectories.mat", 'combined')
