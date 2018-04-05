%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%

function [D_preprocess,D_keep_moving_data] = preProcessing(D,days)

%% delete the not accurate ones and errors
error_index1 = intersect(find(D.gps_accuracy>1000),find(D.missing_data==0)); % Keep the missing data!!!!
D(error_index1,:)=[];
fprintf('delete not accurate data: %i\n',size(error_index1,1));
error_index2 = intersect(find(D.gps_accuracy==0),find(D.missing_data==0));
D(error_index2,:)=[];
fprintf('delete error gps data: %i\n',size(error_index2,1));

%% delete the not changed locations
diff_latitude = D.gps_latitude(1:end-1)-D.gps_latitude(2:end);
diff_longitude = D.gps_longitude(1:end-1)-D.gps_longitude(2:end);
dup_latitude = (roundn(diff_latitude,-6)==0);
dup_longitude = (roundn(diff_longitude,-6)==0);
non_missing = (D.missing_data(2:end)==0); % keep the missing data
dup_index = find(dup_latitude & dup_longitude & non_missing)+1;
D(dup_index,:)=[];
fprintf('delete not changed data: %i\n',size(dup_index,1));

%% add initial orders
order=1:size(D,1);
order=order';
D = [D,table(order)];
D_keep_moving_data=D;

%% seperate days
length_time=size(D,1);
[~,~,d] = ymd(D.time(1));
day_labels = ones(length_time,1);
current_day=1;
for i=2:length_time
   [~,~,d1] = ymd(D.time(i));
   if d1~=d
      d=d1;
      current_day=current_day+1;
   end
   day_labels(i) = current_day;
end
D = [D,table(day_labels)];
fprintf('add day_labels to D\n');
% figure
% index = 1:length_time;
% scatter(index',clusters_final,5,day_labels');
% title('location data of different date in colors');
% fprintf('seperat days: %i\n',current_day);
if days<current_day 
    index_outDate = find(day_labels>days);
    D(index_outDate,:)= [];
else
    days = current_day;
end
fprintf('keep days: %i out of %i \n',max(D.day_labels),current_day);

%% calc duration ( durations of un-changed data are considered)
time_leaving = D.time(2:end);
time_leaving(end+1) = time_leaving(end);
duration = minutes(D.time(2:end)-D.time(1:end-1));
duration(size(D,1))=0;
duration(D.missing_data==1)=0; %% set duration to 0 with missing data
D = [D,table(time_leaving,duration)];
fprintf('calculate leaving time and calculation \n');
%% remove the "on the way" locations
minimal_duaration_minutes =2;
duration_index = find(D.duration<minimal_duaration_minutes & D.duration~=0); % keep missing data
D(duration_index,:)=[];
fprintf('delete small duration locations: %i\n',size(duration_index,1));

%%
D_preprocess = D;
fprintf('only keep data: %i\n',size(D,1));



