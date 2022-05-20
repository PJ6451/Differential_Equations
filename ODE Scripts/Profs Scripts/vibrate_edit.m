% Generating solutions to wave eq for a string
% attached at the ends

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',14);

format compact;
clear
clc
fs=16;
fs2=14;

L = 10;				% length of string
tfin = 20;			% final time
v = 1;				% c value

NptsX=101;
NptsT=101;
Nmodes=6;
Nskip = 2;
x=linspace(0,L,NptsX);
t=linspace(0,tfin,NptsT);
[X,T]=meshgrid(x,t);

a=zeros(1,Nmodes);
b=zeros(1,Nmodes);

for n=1:Nmodes
    int   = -L^2*(n*pi*cos(n*pi)-sin(n*pi))/((pi^2)*(n^2));
    a(n) = (2/L)*int;
    b(n) = (2/(n*pi*v))*int;
end
i=1;
for n=1:Nskip:Nmodes
    ii = n;
    Un = sin(ii*pi/(L)*X).*(a(n).*sin((ii*pi/(L))*v*T) ...
        + b(n).*cos((ii*pi/(L))*v*T));

    figure(1)
    set(gca,'FontSize',fs2);
    subplot(3,1,i)
    set(gca,'FontSize',fs2);
    plot(x./L,Un);
    lab = ['$u_{',num2str(n),'}(x,t)$'];
    if(n==Nmodes)
       xlabel('$x$','Fontsize',fs2); 
    end
    ylabel(lab,'Fontsize',fs2); 
    set(gca,'FontSize',fs2);
    i=i+1;
end


L = 5;				% length of string
tfin = 20;			% final time
v = 1;				% c value

NptsX=101;
NptsT=101;
Nmodes=3;

x=linspace(0,L,NptsX);
t=linspace(0,tfin,NptsT);
[X,T]=meshgrid(x,t);

a=zeros(1,Nmodes);
b=zeros(1,Nmodes);

for n=1:Nmodes
    int  = -4*L^2*(pi*(n - 1/2)*cos(n*pi - pi/2) - sin(n*pi - pi/2))/pi^2/(2*n - 1)^2;
    a(n) = (2/L)*int;
    b(n) = (2/((n-1/2)*pi*v))*int;
end


for n=1:Nmodes
    ii = 2*n-1;
    Un = sin(ii*pi/(2*L)*X).*(a(n).*sin((ii*pi/(2*L))*v*T) ...
        + b(n).*cos((ii*pi/(2*L))*v*T));

    figure(2)
    set(gca,'FontSize',fs2);
    subplot(3,2,2*n-1)
    set(gca,'FontSize',fs2);
    plot(x./L,Un);
    lab = ['$u_{',num2str(n),'}(x,t)$'];
    if(n==Nmodes)
       xlabel('$x$','Fontsize',fs2); 
    end
    ylabel(lab,'Fontsize',fs2); 
    set(gca,'FontSize',fs2);
    
    subplot(3,2,n*2)
    surf(X./L,v.*T./L,Un);
%     colormap(spring)
    xlabel('$x$','Fontsize',fs); ylabel('$t$','Fontsize',fs); zlabel('$u(x,t)$','Fontsize',fs);axis tight
    set(gca,'FontSize',[fs]);
end