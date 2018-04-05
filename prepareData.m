%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%

function D = prepareData(file_addr)

gpsU = importfileGPS(file_addr);
length = size(gpsU,1);
timeUnix = gpsU{2:length,1};
time = datetime(timeUnix, 'ConvertFrom', 'posixtime','TimeZone','America/New_York');
gps_latitude = gpsU{2:length,3};
gps_longitude = gpsU{2:length,4};
gps_accuracy = gpsU{2:length,2};

%% check missing data
% for StudentLife dataset: 20 mins or 10 mins per point
data_original_size = size(time,1);
time_slot =  minutes(time(2:data_original_size)-time(1:data_original_size-1));
missing_data_index = time_slot>30;
missing_data = zeros(data_original_size,1);
missing_data(missing_data_index)=1;
D_original = table(time,gps_latitude,gps_longitude,gps_accuracy,missing_data);
D=D_original;
fprintf('get data: %i\n' ,length);