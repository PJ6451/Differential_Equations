function nyquist_neuralnet

close all;
clear;
clc;

tau = 1/sqrt(3); %somewhere between 1.59157-1.59158 the crossover happens
q1 = 2;
p1 = [1 2 1];
L1 = tf(q1,p1,'InputDelay',tau);

nyquist(L1);

set(gca,'FontSize',12);

