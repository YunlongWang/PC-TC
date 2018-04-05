%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%
function [popularPlacesScore,distInterCluster,popularPlacesClusterIndex,D_cluster]=calculateClusterScore2(distTree,D,max_cluster,durationPerHit_threshold,hitRate_threshold,dayNum)

clusters_final=cluster(distTree,'maxclust',max_cluster);

%% delete the un-changed state
length_clusters_final = size(clusters_final,1);
diff_clusters_final = clusters_final(1:length_clusters_final-1)-clusters_final(2:length_clusters_final);
dup_clusters = (roundn(diff_clusters_final,-4)==0);
dup_index = find(dup_clusters)+1;
D_cluster = [D,table(clusters_final)];
clusters_final(dup_index) = [];
% fprintf('delete not changed clusters: %i\n',size(dup_index,1));

%% calculate durations in each cluster
clusterNum = max(clusters_final);
durationEach = zeros(clusterNum,1);
hitEach=zeros(clusterNum,1);
hitDays=zeros(clusterNum,1);
centerLatEach=zeros(clusterNum,1);
centerLonEach=zeros(clusterNum,1);

for i=1:max(clusters_final)
   tempIndex = find(D_cluster.clusters_final==i);
   durationEach(i) = sum(D_cluster.duration(tempIndex));
   hitEach(i)=size(find(clusters_final==i),1);% here we use clusters removed duplicated
   dayLabels = D_cluster.day_labels(tempIndex);
   hitDays(i) = size(unique(dayLabels),1);
   centerLatEach(i)= mean(D_cluster.gps_latitude(tempIndex));
   centerLonEach(i)=mean(D_cluster.gps_longitude(tempIndex));
end
durationPerHitDay = durationEach./double(hitDays);% shall I remove the days only containing missing data?
distInterCluster = 0;
%% find POIs
popularPlacesClusterIndex = intersect(find(durationPerHitDay > durationPerHit_threshold),find(double(hitDays)/double(dayNum) >hitRate_threshold));
popularPlacesNum = size(popularPlacesClusterIndex,1);
popularPlacesScore=popularPlacesNum;

fprintf('POI score: %i\n',popularPlacesNum);
