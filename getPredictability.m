%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%
function predictability = getPredictability(data_with_cluster)

dataFeed = data_with_cluster;
clusters_final = dataFeed.clusters_final;
length_clusters_final = size(clusters_final,1);
diff_clusters_final = clusters_final(1:length_clusters_final-1)-clusters_final(2:length_clusters_final);
dup_clusters = (roundn(diff_clusters_final,-4)==0);
dup_index = find(dup_clusters)+1;
dataFeed(dup_index,:)=[];

if size(unique(dataFeed.clusters_final),1)>2

entropy = lzentropy(dataFeed.clusters_final);
numPOI = length(unique(dataFeed.clusters_final));
syms x
eqn = -x*log2(x)-(1-x)*log2(1-x)+(1-x)*log2(numPOI-1) == entropy;
predictability = real(solve(eqn,x));

else 
    predictability = 0;
end
