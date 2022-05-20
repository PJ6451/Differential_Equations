% Generating solutions to wave eq for a string
% attached at the ends

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',14);

format compact;
L = 5;				% length of string
H = 1;				% Hight of tent IC
tfin = 20;			% final time
v = 1;				% velocity of propagation in the medium

NptsX=101;
NptsT=101;
Nmodes=12;
Nrows=3;			% number of columns
Nskip=2;			% use every Nskip coefs
Ncols=(Nmodes/Nrows)/Nskip;
x=linspace(0,2*L,NptsX);
t=linspace(0,tfin,NptsT);
[X,T]=meshgrid(x,t);


fs=16;
fs2=14;
plotsols=0;
plotfilm=1;
plotmodes=1;

m = H/(L/2);				% slope of tent IC
U0t = (x<L/2).*x*m+(x>=L/2).*(2*H-x*m);	% initial theoretical steady state distribution
b=zeros(1,Nmodes);

% Fourier coefs for T(x,y=0)=tent temperature (see exercise 7.9.23)
for jj=1:Nskip:Nmodes
 if(jj/2==floor(jj/2))
  bn = 0;
 else    
  bn = (8*H/pi^2)*(-1)^((jj-1)/2)/(jj^2);
 end
 b(jj) = bn;
end

% b(1:Nskip:Nmodes)

U=zeros(NptsT,NptsX);
maxUn=[];
delta=1e-90;

for jj=1:Nskip:Nmodes		% build Fourier series
  Un=b(jj)*sin(((2*jj-1)*pi/(2*L))*X).*cos(((2*jj-1)*pi/(2*L))*v*T);
  maxUn=[maxUn;max(max(Un))]; 
  U=U+Un;
  if(plotsols==1)
   figure(3)
   subplot(2*Nrows,Ncols,jj)
   set(gca,'FontSize',[fs]);
   surf(X,T,Un);
   shading interp
   %colormap(gray)
   zl = ['$u_{',num2str(jj),'}(x,t)$'];
   xlabel('$x$','Fontsize',fs); ylabel('$t$','Fontsize',fs); zlabel(zl,'Fontsize',fs);axis tight
   view([12 46])
   if(jj==1)
     ax=axis;
     ax(5)=-ax(6)/3; 
   end
   axis(ax);
   set(gca,'FontSize',[fs]);
  
   subplot(2*Nrows,Ncols,jj+1)
   set(gca,'FontSize',[fs]);
   surf(X,T,Un);
   shading interp
   colormap(jet)
   view([0 90])
   xlabel('$x$','Fontsize',fs); ylabel('$t$','Fontsize',fs); zlabel('$u(x,t)$','Fontsize',fs);axis tight
   set(gca,'FontSize',[fs]);
   %fprintf('Press any key to continue...\n');
   %pause
   end
end

% Now plot the Fourier combo:
% figure(4)
% clf

  subplot(2,1,1)
  set(gca,'FontSize',[fs]);
  surf(X,T,U);
  shading interp
  %colormap(gray)
  xlabel('$x$','Fontsize',fs); ylabel('$t$','Fontsize',fs); zlabel('$u(x,t)$','Fontsize',fs);axis tight
  view([60 32])
  set(gca,'FontSize',[fs]);

  subplot(2,1,2)
  set(gca,'FontSize',[fs]);
  surf(X,T,U);
  shading interp
  colormap(jet)
  view([0 90])
  xlabel('$x$','Fontsize',fs); ylabel('$t$','Fontsize',fs); zlabel('$u(x,t)$','Fontsize',fs);axis tight
  colorbar
  set(gca,'FontSize',[fs]);

if(plotfilm==1)
 figure(5)
 clf
 xl=min(min(x));
 xr=max(max(x));
 yt=ceil(max(max(U)));
 yb=min(min(U));
 ax=[xl xr yb yt];
 for jj=1:length(U(:,1))
  figure(5)
  set(gca,'FontSize',[fs]);
  plot(x,U(jj,:),'LineWidth',3)
  axis(ax);
  xlabel('$x$','Fontsize',fs);ylabel('$u(x,t)$','Fontsize',fs); 
  mytext=['t=',num2str(t(jj)),'/',num2str(tfin)];
  text(xl+0.1*(xr-xl),yt-0.1*(yt-yb),mytext,'Fontsize',fs)
  
  if(plotmodes==1)
   figure(6)
%   clf
   set(gca,'FontSize',[fs2]);
   for imode=1:Nskip:Nmodes
    Un=b(imode)*sin(((imode-0.5)*pi/L)*x).*cos(((imode-0.5)*pi/L)*v*t(jj));
    subplot(Nrows,Ncols,(imode+1)/Nskip)
    hold on
    set(gca,'FontSize',[fs2]);
    plot(x,Un);
    lab=['$u_{',num2str(imode),'}(x,t)$'];
    if(imode==Nmodes)
       xlabel('$x$','Fontsize',fs2); 
    end
    ylabel(lab,'Fontsize',fs2); 
    ax2=[0 L -maxUn((imode+1)/Nskip)-delta maxUn((imode+1)/Nskip)+delta];
    axis(ax2);
    set(gca,'FontSize',[fs2]);
   end
  end
  
  pause(0.001)
  if(jj==1)
   figure(5)  
   set(gca,'FontSize',[fs]);
   hold on
   plot(x,U0t,'r','LineWidth',3)
   fprintf('Position windows and press any key to continue\n')
   hold off
   pause
  end
  drawnow
 end
end