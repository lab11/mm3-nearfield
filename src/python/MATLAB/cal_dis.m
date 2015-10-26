
function [header_distance, distance, unit_distance, unit_header_distance] = cal_dis(pos)

unit = 1;
unit_distance = 1;
distance = 1;
clear sum,distance,unit,unit_distance;

sum = 0;
distance = 0;
for n = 1: 16 - 1
    
    distance(n) = pos(n+1) - pos(n);
    sum = sum + distance(n); 
end
unit = floor(distance(7)/40);
%unit = 7.5;
%unit = distance(1)/23;
for n = 1: 16 - 1
    header_distance(n) = pos(n+1) - pos(n); 
end
header_unit = header_distance(1)/46;
for n = 1: 16 - 1
    unit_header_distance(n) = header_distance(n)/header_unit; 
end

for n = 1: length(pos) - 1
    distance(n) = pos(n+1) - pos(n); 
end
for n = 1: length(pos) - 1
    unit_distance(n) = distance(n)/unit;
end
unit = unit * 40;
end
