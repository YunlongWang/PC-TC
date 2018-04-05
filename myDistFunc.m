%
% Copyright (c) 2018, Yunlong Wang
% All rights reserved. Please read the "license.txt" for license terms.
%
% This code is part of PC-TC(v0.1)
% 
% Contact Info: yunlong.wang.stud@outlook.com
%
function dist = myDistFunc(XI,XJ)

len = size(XJ,1);
dist = zeros(1,len);
for i=1:len
    dist(i) = 40075017 * distance(XI(1),XI(2),XJ(i,1),XJ(i,2)) / 360;
end
