function a = reviewp8
clear
clc
a    = zeros(10,1);
a(1) = 9;
a(2) = 5;

for i = 3:10
    a(i) = a(i-2)*((i-2)^2 - (i-2) + 2)/(5*(i-1)*(i));
end