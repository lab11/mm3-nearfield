clear all
pos(1) = 6580057;
pos(2) = 6885377;
pos(3) = 7099338;
pos(4) = 7365583;
pos(5) = 7697010;
pos(6) = 8080702;
pos(7) = 8438189;
pos(8) = 8756569;
pos(9) = 9048881;
pos(10) = 9419558;
pos(11) = 9646574;
pos(12) = 9925718;
pos(13) = 10322589;
pos(14) = 10562612;
pos(15) = 10907073;
pos(16) = 11160294;
pos(17) = 11374252;
pos(18) = 11588220;
pos(19) = 11802200;
pos(20) = 12016174;
pos(21) = 12230152;
sum = 0;
distance = 0;
for n = 1: 16 - 1
    distance(n) = pos(n+1) - pos(n);
    sum = sum + distance(n); 
end
unit = sum/345;
for n = 1: length(pos) - 1
    distance(n) = pos(n+1) - pos(n); 
end
for n = 1: length(pos) - 1
    unit_distance(n) = distance(n)/unit;
end