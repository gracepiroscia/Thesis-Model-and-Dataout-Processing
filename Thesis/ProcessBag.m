% Author: Grace Piroscia
%
% The following script functions to convert the USYD campus dataset .ros
% bags into matlab readable structs. The relevant topics and readings that
% are stored include:
%   - '/vn100/imu': Linear X,Y,Z acceleration, TimeStamp information
%   (cumulative time passed is also calculated based off the TimeStamp data
%   and stored)
%
%   - '/zio/odometry/rear': Linear X, Y, Z velocity, TimeStamp information
%
%   - '/gmsl/A0/frame_info': Camera TimeStamp information for the purpose of
%   correlating video footage to sensor readings. 

clear;
clc;

%% Load relevant data and extract relevant topics
filePath = 'bags/2018-04-19-15-05-55_Dataset_Year.bag'; 
bag = rosbag(filePath);

imuBag = select(bag, 'Topic', '/vn100/imu');
odomRearBag = select(bag, 'Topic', '/zio/odometry/rear');
A0cameraFrameBag = select(bag, 'Topic', '/gmsl/A0/frame_info');

imuStructs = readMessages(imuBag,'DataFormat','struct');
odomStructs =  readMessages(odomRearBag,'DataFormat','struct');
A0CamStructs = readMessages(A0cameraFrameBag,'DataFormat','struct');

%% Convert to more readable struct with relevant data extracted
firstEpochTime = 0;
for i = 1:length(imuStructs)
   
    if i == 1
        firstEpochTime = imuStructs{i,1}.Header.Stamp.Sec;
    end
    
    SecondsAccumulated = double(imuStructs{i,1}.Header.Stamp.Sec - firstEpochTime);
    NanoSecondAccumulated  = double(imuStructs{i,1}.Header.Stamp.Nsec);
    
    s.imu(i).TimeSinceEpoch_sec = imuStructs{i,1}.Header.Stamp.Sec;
    s.imu(i).NsecPassedCurrentSec = imuStructs{i,1}.Header.Stamp.Nsec;
    s.imu(i).Time_since_start_secs = SecondsAccumulated + NanoSecondAccumulated*1e-9;
    s.imu(i).LinearAcc_X = imuStructs{i,1}.LinearAcceleration.X;
    s.imu(i).LinearAcc_Y = imuStructs{i,1}.LinearAcceleration.Y;
    s.imu(i).LinearAcc_Z = imuStructs{i,1}.LinearAcceleration.Z;
    
end

firstEpochTime = 0;
for i = 1:length(odomStructs)
   
    if i == 1
        firstEpochTime = odomStructs{i,1}.Header.Stamp.Sec;
    end
    
    SecondsAccumulated = double(odomStructs{i,1}.Header.Stamp.Sec - firstEpochTime);
    NanoSecondAccumulated  = double(odomStructs{i,1}.Header.Stamp.Nsec);
    
    s.odomRear(i).TimeSinceEpoch_sec = odomStructs{i,1}.Header.Stamp.Sec;
    s.odomRear(i).NsecPassedCurrentSec = odomStructs{i,1}.Header.Stamp.Nsec;
    s.odomRear(i).Time_since_start_secs = SecondsAccumulated + NanoSecondAccumulated*1e-9;
    s.odomRear(i).LinearVel_X = odomStructs{i,1}.Twist.Twist.Linear.X;
    s.odomRear(i).LinearVel_Y = odomStructs{i,1}.Twist.Twist.Linear.Y;
    s.odomRear(i).LinearVel_Z = odomStructs{i,1}.Twist.Twist.Linear.Z;
    
end

firstEpochTime = 0;
for i = 1:length(A0CamStructs)
   
    if i == 1
        firstEpochTime = A0CamStructs{i,1}.Header.Stamp.Sec;
    end
    
    SecondsAccumulated = double(A0CamStructs{i,1}.Header.Stamp.Sec - firstEpochTime);
    NanoSecondAccumulated  = double(A0CamStructs{i,1}.Header.Stamp.Nsec);
    
    s.Camera(i).TimeSinceEpoch_sec = A0CamStructs{i,1}.Header.Stamp.Sec;
    s.Camera(i).NsecPassedCurrentSec = A0CamStructs{i,1}.Header.Stamp.Nsec;
    s.Camera(i).Time_since_start_secs = SecondsAccumulated + NanoSecondAccumulated*1e-9;
    
end

%% Save struct
FileName = filePath(6:end-4); %removing .bag extension and parent directory bags/
FileName = FileName + ".mat";
save(FileName, 's')