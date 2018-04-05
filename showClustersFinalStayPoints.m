%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%

function showClustersFinalStayPoints(clusters_final,data_preprecessing_stayPoints,title_str)

cluster_final_DBSCAN = zeros(size(data_preprecessing_stayPoints,1),1);
for i=1:max(clusters_final)
   index_1 = find(clusters_final==i);
   index_2 = ismember(data_preprecessing_stayPoints.stayPoints_labels,index_1);
   cluster_final_DBSCAN(index_2)=i;
end

poiPosition = find(cluster_final_DBSCAN~=0);
figure,       
[~,colorCode] = ismember(cluster_final_DBSCAN,unique(cluster_final_DBSCAN));       
scatter(data_preprecessing_stayPoints.gps_longitude(poiPosition),data_preprecessing_stayPoints.gps_latitude(poiPosition),30,colorCode(poiPosition),'filled');
title([title_str num2str(max(clusters_final))]);
plot_google_map;
