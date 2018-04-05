%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%
function [D_cluster_init,POIsIndex,best_cluster_size,realPOIIndex] = findBestCluster(D,localPOI,hitRate_local)

%% clustering by different cutoff for daily pattern
pairs = [D.gps_latitude D.gps_longitude];

cluster_min=2;
cluster_max=size(unique(D.day_labels),1)*2;

iterate = 1;
durationPerHitDay_threshold=120; % minutes for global POI
hitRate_threshold=0.63; % for global POI

if localPOI ==1 % for local POI
    cluster_max=size(D,1);
    durationPerHitDay_threshold = 30;
    hitRate_threshold = hitRate_local; 
end

D_old = D;

%coordinate mapping
for i=1:size(pairs,1)
    pairs(i,2) = pairs(i,2)*cos( degtorad(pairs(i,1)));
end

%% hierarchical clustering
distTree = linkage(pairs,'complete');

window=1;
popularPlacesScore = zeros(floor((cluster_max-cluster_min+1)/window),1);
distInterCluster = zeros(floor((cluster_max-cluster_min+1)/window),1);

cluster_history = zeros(floor((cluster_max-cluster_min+1)/window),1);
max_cluster = cluster_min;
best_cluster = 1;
best_score_temp = 0;
while max_cluster <= cluster_max
    D = D_old;
    fprintf('-----\nclustering number: %i\n',max_cluster);
    [popularPlacesScore(iterate),distInterCluster(iterate),~,~]= calculateClusterScore2(distTree,D,max_cluster,durationPerHitDay_threshold,hitRate_threshold,max(D.day_labels));
    if popularPlacesScore(iterate)>best_score_temp
        best_score_temp = popularPlacesScore(iterate);
        best_cluster = iterate;
    end
    keep_lower_state=false;
    if best_score_temp >= popularPlacesScore(iterate)
        if length(find(popularPlacesScore(best_cluster:iterate)<=best_score_temp))>50
            keep_lower_state=true;
        end
    end

    cluster_history(iterate)=max_cluster;
    
    % stop when detect these states
    if popularPlacesScore(iterate)==0 ||  best_score_temp-popularPlacesScore(iterate)> 10 || keep_lower_state==true
       break; 
    end
    
    iterate=iterate+1;
    max_cluster = max_cluster+window;

end

%% using the best cluster
[best_score,~] = max(popularPlacesScore);
best_cluster_size=0;
if best_score ~=0
    best_score_index_all = find(popularPlacesScore==best_score);
    best_score_index = best_score_index_all(end);
    best_cluster_size=cluster_history(best_score_index);
end

D = D_old;
fprintf('-----\n best clustering number: %i\n',best_cluster_size);

if best_cluster_size ==0
    POIsIndex=1;
    clusters_final = ones(size(D,1),1);
    realPOIIndex = clusters_final;
    D_cluster_init = [D,table(clusters_final)];
else
    [~,~,realPOIIndex,~]= calculateClusterScore2(distTree,D,best_cluster_size,durationPerHitDay_threshold,hitRate_threshold,max(D.day_labels));

    durationPerHitDay_threshold=30; % to include all POIs
    hitRate_threshold = 0.13;

    [~,~,POIsIndex,D_cluster_init]= calculateClusterScore2(distTree,D,best_cluster_size,durationPerHitDay_threshold,hitRate_threshold,max(D.day_labels));
end
