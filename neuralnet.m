function neuralnet

close all;
clear;
clc;

% Parameters
a_1 = -2;
a_2 = 1;
t1 = pi/3.5;
t2 = pi/3.5;
% t1 = .1;
% t2 = .2;

history = [20.0;15.0];
tspan   = [0 100];
opts    = ddeset('RelTol',1e-5,'AbsTol',1e-8);
% Solve the ODEs that arise when there is no delay.
sol1 = dde23(@ddefun,[],history,tspan,[],a_1,a_2);

% Solve the DDEs that arise when there is a delay of r.
sol2 = dde23(@ddefun,[t1;t2],history,tspan,opts,a_1,a_2);

figure(1)
plot(sol1.x,sol1.y(1,:),'b-','Linewidth', 3), hold on;
plot(sol1.x,sol1.y(2,:),'r-','Linewidth', 3), hold on;
plot(sol2.x,sol2.y(1,:),'g-','Linewidth', 3), hold on;
plot(sol2.x,sol2.y(2,:),'k-','Linewidth', 3), hold on;
xlabel('t');
ylabel('x(t)');
legend({'No Delay','No Delay','Delay','Delay'},'Location','NorthEast');
set(gca,'FontSize',12);
grid on;


% equation being solved
function dudt = ddefun(t,x,Z,a_1,a_2) 
if isempty(Z)     % ODEs
   dudt = [-x(1) + a_1*tanh(x(2)); -x(2) + a_2*tanh(x(1))];
else
  xlag1 = Z(:,1);
  xlag2 = Z(:,2);
  dudt = [-x(1) + a_1*tanh(xlag2(2)); -x(2) + a_2*tanh(xlag1(1))];
end
