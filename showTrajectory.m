%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%

function showTrajectory(D,note)

ifShow = true;
if ifShow == true
    
    % modified jet-colormap
    temp_data = D.gps_latitude;
    cd_3 = colormap('jet'); % take your pick (doc colormap)
    cd_3 = interp1(linspace(min(temp_data),max(temp_data),length(cd_3)),cd_3,temp_data); % map color to one dimention value
    cd_3 = uint8(cd_3'*255); % need a 4xN uint8 array
    cd_3(4,:) = 255; % last column is transparency
    
    p1=geoshow(D.gps_latitude,D.gps_longitude,'DisplayType','point','Marker','.','Color','black','MarkerEdgeColor','auto','MarkerSize',10);
    hold on
    p2=geoshow(D.gps_latitude,D.gps_longitude,'DisplayType','line','Color','red');
    title(note);
    xlabel('latitude');ylabel('longitude');
    
    drawnow;
    set(p1.Edge, 'ColorBinding','interpolated', 'ColorData',cd_3);
    set(p2.Edge, 'ColorBinding','interpolated', 'ColorData',cd_3);
    plot_google_map;
end
