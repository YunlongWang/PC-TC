%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is to demonstrate the PC-TC algorithm in comparison with other
% approaches.
%
% Project Name: PC-TC(v1.0)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%
close all;
clc;
clear;

%% init variables
userSize = 1;
durationPercentage = zeros(userSize,5);
POI_number = zeros(userSize,5);
predictability= zeros(userSize,5);
computation_time = zeros(userSize,5);
stay_point_detetion_time = zeros(userSize,1);
data_size_preprocessing = zeros(userSize,1);

best_cluster_num = zeros(userSize,1);
best_cluster_num_local = zeros(userSize,1);
real_global_POI_num = zeros(userSize,1);
predictability_global = zeros(userSize,1);
missing_duration = zeros(userSize,1);

%% parameters
distThretSP= 50;% meters
timeThretSP = 30;% minutes
hitRate_local = 0.6; % F_vd, change this to 0.1-0.9

ifShowFigures = true;

%% import and prepare data
u=1;
% the data is from StudentLife dataset (http://studentlife.cs.dartmouth.edu/)
data_original = prepareData('gps_u16.csv');

%% pre-processing
daysDuration = 60;
[data_preprocessing,data_keep_moving] = preProcessing(data_original,daysDuration);
total_time = minutes(data_preprocessing.time(end)-data_preprocessing.time(1));
data_size_preprocessing(u,1)=size(data_preprocessing,1);

missing_data_index = find(data_preprocessing.missing_data==1);
if missing_data_index(end)==length(data_preprocessing.missing_data)
    missing_data_index(end)=[];
end
missing_data_original_index = data_preprocessing.order(missing_data_index);
missing_duration(u) = sum(minutes(data_original.time(missing_data_original_index+1)-data_original.time(missing_data_original_index)));

%% data preview in 2D and 3D
if ifShowFigures==true
    showTrajectory(data_preprocessing,'');
    figure,
    xx = data_preprocessing.gps_longitude;
    yy = data_preprocessing.gps_latitude;
    zz = data_preprocessing.time;
    p=plot3(xx,yy,zz);
    xlabel('latitude');ylabel('longitude');zlabel('time');
    % modified jet-colormap
    nn = size(data_preprocessing,1);
    cd = [uint8(jet(nn)*255) uint8(ones(nn,1))].';
    temp_data = yy;
    cd_3 = colormap('jet'); % take your pick (doc colormap)
    cd_3 = interp1(linspace(min(temp_data),max(temp_data),length(cd_3)),cd_3,temp_data); % map color to one dimention value
    cd_3 = uint8(cd_3'*255); % need a 4xN uint8 array
    cd_3(4,:) = 255; % last column is transparency
    drawnow;
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd);
end

%% get stay points

% coordinate mapping
stayTic = tic;
[stayPoints,data_preprecessing_stayPoints]= getStayPoints(data_preprocessing,distThretSP,timeThretSP);
data_stay = [stayPoints.lat stayPoints.lng];

for i=1:size(data_stay,1)
    data_stay(i,2) = data_stay(i,2)*cos( degtorad(data_stay(i,1)));
end

stay_point_detetion_time(u) = toc(stayTic);

%% OPTICS
optics_tic = tic;
[RD,CD,order] = optics_JClust(data_stay,round(log10(size(data_stay,1))));
%figure,bar(RD)
clusters_OPTICS = extract_clusters_optics(RD,order,round(log10(size(data_stay,1))));
clusters = zeros(size(data_stay,1),1);
POI_number(u,1) = size(clusters_OPTICS,2);
for i=1:size(clusters_OPTICS,2)
    clusters(clusters_OPTICS{1,i})=i;
end
computation_time(u,1)=toc(optics_tic)+stay_point_detetion_time;
clusters_final = clusters;
data_with_clusters = [stayPoints table(clusters_final)];
POIs = formatPOIsWithStayPoints(data_preprocessing,stayPoints,data_with_clusters,clusters);
durationPercentage(u,1) = sum(POIs.POIs_average_duration .* POIs.POIs_visit_times)/(total_time-missing_duration(u));
predictability(u,1)=getPredictability(data_with_clusters);

if ifShowFigures==true
    showClustersFinalStayPoints(data_with_clusters.clusters_final,data_preprecessing_stayPoints,'OPTICS');
end

%% DBSCAN
dbscan_tic = tic;
[D_with_cluster_SOC,POIs_SOC]=findPOIs_DBSCAN(stayPoints,data_preprocessing,distThretSP);
computation_time(u,2)=toc(dbscan_tic)+stay_point_detetion_time;
POI_number(u,2) = size(POIs_SOC,1);
durationPercentage(u,2) = sum(POIs_SOC.POIs_average_duration .* POIs_SOC.POIs_visit_times)/(total_time-missing_duration(u));
predictability(u,2)=getPredictability(D_with_cluster_SOC);

if ifShowFigures==true
    showClustersFinalStayPoints(D_with_cluster_SOC.clusters_final,data_preprecessing_stayPoints,'DBSCAN');
end

%% DB
myfunc = @(X,K)(clusterdata(X,'maxclust',K,'linkage','complete'));
db_tic = tic;
eva = evalclusters(data_stay,myfunc,'DaviesBouldin','KList',1:floor(size(data_stay,1)/2));
%figure,plot(eva);
clusters = clusterdata(data_stay,'maxclust',eva.OptimalK);
computation_time(u,3)=toc(db_tic)+stay_point_detetion_time;
clusters_final = clusters;
POI_number(u,3)= eva.OptimalK;
data_with_clusters = [stayPoints table(clusters_final)];
POIs = formatPOIsWithStayPoints(data_preprocessing,stayPoints,data_with_clusters,clusters);
durationPercentage(u,3) = sum(POIs.POIs_average_duration .* POIs.POIs_visit_times)/(total_time-missing_duration(u));
predictability(u,3)=getPredictability(data_with_clusters);

if ifShowFigures==true
    showClustersFinalStayPoints(data_with_clusters.clusters_final,data_preprecessing_stayPoints,'DB');
end

%% SC
sc_tic = tic;
eva = evalclusters(data_stay,myfunc,'silhouette','KList',1:floor(size(data_stay,1)/2));
%figure,plot(eva);
clusters = clusterdata(data_stay,'maxclust',eva.OptimalK);
computation_time(u,4)=toc(sc_tic)+stay_point_detetion_time;
clusters_final = clusters;
POI_number(u,4)= eva.OptimalK;
data_with_clusters = [stayPoints table(clusters_final)];
POIs = formatPOIsWithStayPoints(data_preprocessing,stayPoints,data_with_clusters,clusters);
durationPercentage(u,4) = sum(POIs.POIs_average_duration .* POIs.POIs_visit_times)/(total_time-missing_duration(u));
predictability(u,4)=getPredictability(data_with_clusters);

if ifShowFigures==true
    showClustersFinalStayPoints(data_with_clusters.clusters_final,data_preprecessing_stayPoints,'SC');
end

%% PC-TC
mytic = tic;
% POIType: 0-global, 1-local
POIType=0;
[D_cluster_init,POI_cluster,best_cluster_num(u),realPOIIndex] = findBestCluster(data_preprocessing,POIType,hitRate_local);
real_global_POI_num(u)=size(unique(realPOIIndex),1);
totalClusterSize = max(D_cluster_init.clusters_final);
timeGlobalPOI = toc(mytic);

% predictability
calcPredictability = true;
if calcPredictability == true
    dataFeed = D_cluster_init;
    [poiPosition,~] = ismember(dataFeed.clusters_final,POI_cluster);
    dataFeed = dataFeed(poiPosition,:);
    clusters_final = dataFeed.clusters_final;
    length_clusters_final = size(clusters_final,1);
    diff_clusters_final = clusters_final(1:length_clusters_final-1)-clusters_final(2:length_clusters_final);
    dup_clusters = (roundn(diff_clusters_final,-4)==0);
    dup_index = find(dup_clusters)+1;
    dataFeed(dup_index,:)=[];
    
    if length(POI_cluster)>2
        entropy = lzentropy(dataFeed.clusters_final);
        numPOI = length(unique(dataFeed.clusters_final));
        syms x
        eqn = -x*log2(x)-(1-x)*log2(1-x)+(1-x)*log2(numPOI-1) == entropy;
        solx = solve(eqn,x);
        predictability_global(u)= solx;
    end
end

% visualization
if ifShowFigures==true
    figure,
    [poiPosition,~] = ismember(D_cluster_init.clusters_final,POI_cluster);
    [~,colorCode] = ismember(D_cluster_init.clusters_final,unique(D_cluster_init.clusters_final));
    scatter(D_cluster_init.gps_longitude(poiPosition),D_cluster_init.gps_latitude(poiPosition),100,colorCode(poiPosition));
    title(['Global POIs' num2str(size(unique(D_cluster_init.clusters_final(poiPosition)),1))]);
    plot_google_map;
end

% cluster statistics
mytic = tic;
[cluster_seq, GlobalPOIs, GlobalPOIs_visit] = formatPOIs(D_cluster_init,POI_cluster);

% go for local POIs
local_POI_cluster = cell(size(POI_cluster,1),2);
D_after_clustering_local= cell(size(POI_cluster,1),2);
LocalPOIs= cell(size(POI_cluster,1),2);
LocalPOIs_visit= cell(size(POI_cluster,1),2);
localRepeatedTrajectories= cell(size(POI_cluster,1),2);
localRoutinePatterns= cell(size(POI_cluster,1),2);

global_POI_id = GlobalPOIs.popularPlacesClusterIndex;
final_POIs = [GlobalPOIs table(global_POI_id)];
final_data_clusters = D_cluster_init; % iniit fianl_data_clusters
POIType=1;
for i=1:size(POI_cluster,1)
    
    globalPOI = POI_cluster(i);
    
    localData = D_cluster_init(D_cluster_init.clusters_final==globalPOI,:);%
    
    localLatLng = [localData.gps_latitude localData.gps_longitude];
    uniquePos = unique(localLatLng,'rows');
    LocalPOIs{i,2}=globalPOI;
    if size(uniquePos,1)>1
        
        localData.clusters_final=[];
        [D_local_cluster_init,local_POI_cluster{i,1},bestClusterNumLocal,~] = findBestCluster(localData,POIType,hitRate_local);
        if bestClusterNumLocal > best_cluster_num_local(u)
            best_cluster_num_local(u)=bestClusterNumLocal;
        end
        % combine the local cluster into global clusters
        for j=1:size(D_local_cluster_init,1)
            if size(local_POI_cluster{i,1},1)>1 % only if there are more than 1 local POIs, we update
                index = find(final_data_clusters.order==D_local_cluster_init.order(j));
                final_data_clusters.clusters_final(index) = D_local_cluster_init.clusters_final(j)+totalClusterSize;
            end
        end
        
        local_POI_cluster{i,2} = globalPOI;
        local_POI_cluster{i,3} = totalClusterSize; %% poi number base
        
        % more than one cluster
        if size(local_POI_cluster{i,1},1) > 1
            [D_after_clustering_local{i,1}, LocalPOIs{i,1}, LocalPOIs_visit{i,1}] = formatPOIs(D_local_cluster_init,local_POI_cluster{i,1});
            D_after_clustering_local{i,2}=globalPOI;
            LocalPOIs_visit{i,2}=globalPOI;
            % combine local POIs into GloablePOIs
            temp = LocalPOIs{i,1};
            temp.popularPlacesClusterIndex = temp.popularPlacesClusterIndex+ totalClusterSize;
            global_POI_id = ones(size(local_POI_cluster{i,1},1),1)*globalPOI;
            temp = [temp table(global_POI_id)];
            final_POIs = vertcat(final_POIs,temp);
        end
        
        %update totoal cluster number
        totalClusterSize = totalClusterSize+max(D_local_cluster_init.clusters_final);
        
    end
    
end

computation_time(u,5)=toc(mytic)+timeGlobalPOI;


% calculate the time spend percentage in POIs
pois_u_num=0;
pois_u_duration=0;
for i=1:size(LocalPOIs,1)
    if ~isempty(LocalPOIs{i,1})
        pois_u_num = pois_u_num + size(LocalPOIs{i,1}.popularPlacesClusterIndex,1);
        pois_u_duration = pois_u_duration + sum(LocalPOIs{i,1}.POIs_visit_times .* LocalPOIs{i,1}.POIs_average_duration);
    else
        pois_u_num = pois_u_num + 1;
        global_poi_index = find(GlobalPOIs.popularPlacesClusterIndex == LocalPOIs{i,2});
        pois_u_duration = pois_u_duration + GlobalPOIs.POIs_visit_times(global_poi_index) *  GlobalPOIs.POIs_average_duration(global_poi_index);
    end
end
POI_number(u,5)=pois_u_num;
durationPercentage(u,5)=pois_u_duration/(total_time-missing_duration(u));

% predictability
if calcPredictability == true
    dataFeed = final_data_clusters;
    [poiPosition,loc] = ismember(dataFeed.clusters_final,final_POIs.popularPlacesClusterIndex);
    dataFeed = dataFeed(poiPosition,:);
    clusters_final = dataFeed.clusters_final;
    length_clusters_final = size(clusters_final,1);
    diff_clusters_final = clusters_final(1:length_clusters_final-1)-clusters_final(2:length_clusters_final);
    dup_clusters = (roundn(diff_clusters_final,-4)==0);
    dup_index = find(dup_clusters)+1;
    dataFeed(dup_index,:)=[];
    
    if length(final_POIs.popularPlacesClusterIndex)>2
        entropy = lzentropy(dataFeed.clusters_final);
        numPOI = length(unique(dataFeed.clusters_final));
        syms x
        eqn = -x*log2(x)-(1-x)*log2(1-x)+(1-x)*log2(numPOI-1) == entropy;
        solx = solve(eqn,x);
        predictability(u,5)= solx;
    end
end

% visualization
if ifShowFigures ==true
    figure,
    [~,colorCode] = ismember(final_data_clusters.clusters_final,unique(final_data_clusters.clusters_final));
    scatter(final_data_clusters.gps_longitude(poiPosition),final_data_clusters.gps_latitude(poiPosition),30,colorCode(poiPosition),'filled');
    title(['local POIs' num2str(size(unique(final_data_clusters.clusters_final(poiPosition)),1))]);
    plot_google_map;
end

%%
if ifShowFigures ==true
    
    methods_str = {'OPTICS','DBSCAN','DB','SC','PC-TC'};
    methods_global_str={'OPTICS','DBSCAN','DB','SC','PC-TC','Global POIs'};
    fontSize = 8;
    
    figure,bar([POI_number real_global_POI_num],'DisplayName','numbers of POIs');
    title('Numbers of POIs');xlabel('user');
    set(gca,'XTick',1:6,'XTickLabel',methods_global_str,'FontSize', fontSize);   
    
    figure,bar([predictability predictability_global],'DisplayName','predictability limit');
    title('predictability limit');xlabel('user');
    set(gca,'XTick',1:6,'XTickLabel',methods_global_str,'FontSize', fontSize);
    
    figure;
    bar(computation_time,'DisplayName','calculate time');
    title('computation time');xlabel('user');ylabel('seconds');
    set(gca,'XTick',1:6,'XTickLabel',methods_global_str,'FontSize', fontSize);
    
end


