function protein

close all;
clear;
clc;

% Parameters
alpha   = 100;
beta    = 1.1;
gamma   = 1;
tau     = 50;
history = [20.0];
tspan   = [0 600];
opts    = ddeset('RelTol',1e-5,'AbsTol',1e-8);
% Solve the ODEs that arise when there is no delay.
sol1 = dde23(@ddefun,[],history,tspan,[],alpha, beta, gamma);

% Solve the DDEs that arise when there is a delay of r.
sol2 = dde23(@ddefun,[tau],history,tspan,opts,alpha,beta,gamma);

figure(1)
plot(sol1.x,sol1.y,'b-','Linewidth', 3), hold on;
plot(sol2.x,sol2.y,'k-','Linewidth', 3), hold on;
xlabel('t');
ylabel('x(t)');
legend('No Delay','Delay','Location','NorthEast');
set(gca,'FontSize',12);
grid on;


% equation being solved
function dpdt = ddefun(t,x,Z,alpha,beta,gamma) 
if isempty(Z)     % ODEs
   dpdt = alpha - beta*x(1) - gamma*x(1);
else
  xlag = Z(:,1);
  dpdt = alpha - beta*x(1) - gamma*xlag(1);
end

