function nyquist_protein

close all;
clear;
clc;

tau   = 10.0;
gamma = 1;
beta  = 1.1;

q1 = gamma;
p1 = [1 beta];
L1 = tf(q1,p1,'InputDelay',tau);

nyquist(L1);

set(gca,'FontSize',12);

