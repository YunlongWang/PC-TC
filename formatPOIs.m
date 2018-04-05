%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%
function [cluster_seq,POIs,POIs_visit] = formatPOIs(D_cluster_init,popularPlacesClusterIndex)

clusters_final=D_cluster_init.clusters_final;

%% calc durations for each cluster
length_clusters_final = size(clusters_final,1);
diff_clusters_final = clusters_final(1:length_clusters_final-1)-clusters_final(2:length_clusters_final);
cluster_start = find(roundn(diff_clusters_final,-4)~=0)+1;
cluster_start = [1;cluster_start];

dup_clusters = (roundn(diff_clusters_final,-4)==0);
dup_index = find(dup_clusters)+1;
cluster_seq = D_cluster_init;
cluster_seq(dup_index,:) = [];

diff_clusters_final = clusters_final(2:length_clusters_final)-clusters_final(1:length_clusters_final-1);
cluster_end = find(roundn(diff_clusters_final,-4)~=0);
cluster_end = [cluster_end;length_clusters_final];

for i=1:size(cluster_seq,1)
    cluster_seq.time_leaving(i) = D_cluster_init.time_leaving(cluster_end(i));
    cluster_seq.duration(i) = sum(D_cluster_init.duration(cluster_start(i):cluster_end(i)));
end

%% 
clusters_num = max(clusters_final);
durationEach = zeros(clusters_num,1);
hitEach=zeros(clusters_num,1);
latitudeEachCluster = zeros(clusters_num,1);
latitudeRadiusEachCluster = zeros(clusters_num,1);
longitudeEachCluster = zeros(clusters_num,1);
longitudeRadiusEachCluster = zeros(clusters_num,1);
maxRadiusEachCluster = zeros(clusters_num,1);
radiusEachCluster = zeros(clusters_num,1);
maxLatitudeEachCluster = zeros(clusters_num,1);
minLatitudeEachCluster = zeros(clusters_num,1);
maxLongitudeEachCluster = zeros(clusters_num,1);
minLongitudeEachCluster = zeros(clusters_num,1);

clusters_visit = cell(clusters_num,2);

for i=1:max(clusters_final)
   durationEach(i) = sum(D_cluster_init.duration(D_cluster_init.clusters_final==i));
   hitEach(i)=size(find(cluster_seq.clusters_final==i),1); %% here we use cluster_seq, which has removed the duplicated clusters
   allInThisCluster = D_cluster_init(D_cluster_init.clusters_final==i,:);
   latMean= median(D_cluster_init.gps_latitude(D_cluster_init.clusters_final==i));
   latitudeRadiusEachCluster(i)= max(abs(D_cluster_init.gps_latitude(D_cluster_init.clusters_final==i)-latMean));
   latitudeEachCluster(i)=latMean;
   lonMean =  median(D_cluster_init.gps_longitude(D_cluster_init.clusters_final==i));
   longitudeRadiusEachCluster(i)= max(abs(D_cluster_init.gps_longitude(D_cluster_init.clusters_final==i)-lonMean));
   longitudeEachCluster(i) = lonMean;
   maxRadiusEachCluster(i) = max([ latitudeRadiusEachCluster(i)  longitudeRadiusEachCluster(i)]);
   E = referenceEllipsoid('earth');
   indexInDInit = find(D_cluster_init.clusters_final==i);
   dist = zeros(size(indexInDInit,1),1);
   for j=1:size(indexInDInit,1)
       dist(j) = distance(allInThisCluster.gps_latitude(j),allInThisCluster.gps_longitude(j),latMean,lonMean,E);
   end
   radiusEachCluster(i)=max(dist); % in meters
   cluster_visit_time = cluster_seq.time(cluster_seq.clusters_final==i);
   cluster_visit_duration = cluster_seq.duration(cluster_seq.clusters_final==i);
   cluster_visit = table(cluster_visit_time,cluster_visit_duration);
   clusters_visit{i,1}= cluster_visit;
   clusters_visit{i,2}=i;
   
   maxLatitudeEachCluster(i) = max(D_cluster_init.gps_latitude(D_cluster_init.clusters_final==i));
   minLatitudeEachCluster(i) = min(D_cluster_init.gps_latitude(D_cluster_init.clusters_final==i));
   maxLongitudeEachCluster(i) = max(D_cluster_init.gps_longitude(D_cluster_init.clusters_final==i));
   minLongitudeEachCluster(i) = min(D_cluster_init.gps_longitude(D_cluster_init.clusters_final==i));
end

POIs_center_latitude = latitudeEachCluster(popularPlacesClusterIndex);
POIs_center_longitude = longitudeEachCluster(popularPlacesClusterIndex);
POIs_radius = radiusEachCluster(popularPlacesClusterIndex);
POIs_max_latitude =  maxLatitudeEachCluster(popularPlacesClusterIndex);
POIs_min_latitude = minLatitudeEachCluster(popularPlacesClusterIndex);
POIs_max_longitude = maxLongitudeEachCluster(popularPlacesClusterIndex);
POIs_min_longitude = minLongitudeEachCluster(popularPlacesClusterIndex);
POIs_visit_times = hitEach(popularPlacesClusterIndex);
POIs_coverage_days = ones(size(popularPlacesClusterIndex,1),1)*max(D_cluster_init.day_labels);% not D!!! because D may not contain all the day labels
POIs_average_duration = durationEach(popularPlacesClusterIndex)./POIs_visit_times;
POIs_visit = clusters_visit(popularPlacesClusterIndex,:);
POIs = table(popularPlacesClusterIndex,POIs_center_latitude,POIs_center_longitude,POIs_radius,POIs_max_latitude,POIs_min_latitude,POIs_max_longitude,POIs_min_longitude,POIs_visit_times,POIs_coverage_days,POIs_average_duration);
