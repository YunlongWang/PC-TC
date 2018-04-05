%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%

function [stayPoints,data_original_with_stayPoints]=getStayPoints(data_original,distThret,timeThret)

%% find POIs
dataLength = size(data_original,1);
stayPoints_labels=zeros(dataLength,1);
stayPoints_label_num=1;
stayPoints =[];
i=1;
while i<dataLength
   token=0;
   if data_original.missing_data(i)==1 %% avoid starting with a missing data
      i=i+1;
      continue;
   end
   for j=i+1:dataLength
      dist= 40075017 * distance(data_original.gps_latitude(i),data_original.gps_longitude(i),data_original.gps_latitude(j),data_original.gps_longitude(j)) / 360; % meters
      if dist>distThret || data_original.missing_data(j)==1 %% avoid including a missing data
         timeSpan = minutes(data_original.time(j)-data_original.time(i));
         if timeSpan > timeThret
             clusterPoints = data_original(i:j-1,:);
             lat=mean(clusterPoints.gps_latitude);
             lng=mean(clusterPoints.gps_longitude);
             arvT=data_original.time(i);
             levT=data_original.time(j);
             arvN=i;
             levN=j-1;
             day_labels=data_original.day_labels(i);
             missing_data=0;
             if data_original.missing_data(j)==1
                missing_data=1; 
             end
             stayPoints = vertcat(stayPoints,table(lat,lng,arvT,levT,timeSpan,arvN,levN,day_labels,missing_data));
             stayPoints_labels(i:j-1)=stayPoints_label_num;
             stayPoints_label_num = stayPoints_label_num+1;
             i=j;
             token=1;
         end
         break;
      end
   end
   if token==0
       i=i+1;
   end
end

data_original_with_stayPoints = [data_original table(stayPoints_labels)];


