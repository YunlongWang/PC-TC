%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%

function [data_with_clusters_SOC,POIs]=findPOIs_DBSCAN(stayPoints,data_original,epsilon)

data_DBSCAN = [stayPoints.lat stayPoints.lng];
minpts = round(log10(size(stayPoints,1))); % points

%% DBSCAN
IDX=DBSCAN(data_DBSCAN,epsilon,minpts);

%% format data
clusters_final = IDX;
data_with_clusters_SOC = [stayPoints table(clusters_final)];

%% format POIs
maxNum = max(IDX);
latitudeEachCluster = zeros(maxNum,1);
longitudeEachCluster = zeros(maxNum,1);
radiusEachCluster = zeros(maxNum,1);
hitEach=zeros(maxNum,1);
durationEach=zeros(maxNum,1);

for i=1:maxNum
   latitudeEachCluster(i) = mean(stayPoints.lat(IDX==i));
   longitudeEachCluster(i) = mean(stayPoints.lng(IDX==i));
   maxLat = max(stayPoints.lat(IDX==i));
   minLat = min(stayPoints.lat(IDX==i));
   if maxLat-latitudeEachCluster(i)>latitudeEachCluster(i)-minLat
      radiusLat = maxLat-latitudeEachCluster(i);
   else
      radiusLat = latitudeEachCluster(i)-minLat;
   end
   maxLng = max(stayPoints.lng(IDX==i));
   minLng = min(stayPoints.lng(IDX==i));
   if maxLng-longitudeEachCluster(i)>longitudeEachCluster(i)-minLng
      radiusLng = maxLng-longitudeEachCluster(i);
   else
      radiusLng = longitudeEachCluster(i)-minLng;
   end
   
   if radiusLat >radiusLng
       distInDegree = radiusLat;
   else
       distInDegree = radiusLng;
   end
   
   radiusEachCluster(i)=40075017 * distance(latitudeEachCluster(i),longitudeEachCluster(i),latitudeEachCluster(i)+distInDegree,longitudeEachCluster(i))/360;
   
end

for i=1:maxNum
   hitEach(i)=size(find(data_with_clusters_SOC.clusters_final==i),1);
   durationEach(i)=sum(data_with_clusters_SOC.timeSpan(data_with_clusters_SOC.clusters_final==i));
end

popularPlacesClusterIndex = (1:maxNum)';
POIs_center_latitude = latitudeEachCluster;
POIs_center_longitude = longitudeEachCluster;
POIs_radius = radiusEachCluster;
POIs_visit_times = hitEach;
POIs_coverage_days = ones(maxNum,1)*max(data_original.day_labels);% not D!!! because D may not contain all the day labels
POIs_average_duration = durationEach./POIs_visit_times;
POIs = table(popularPlacesClusterIndex,POIs_center_latitude,POIs_center_longitude,POIs_radius,POIs_visit_times,POIs_coverage_days,POIs_average_duration);


