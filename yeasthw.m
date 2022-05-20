clear
clc
load yeast.mat

%part a
Nc = yeast{end,2};
N0 = yeast{1,2};
N  = N0/Nc;
A = zeros(19,1);
for i = 1:size(yeast,1)
    A(i)  = ((1-N)/N)*((yeast{i,2})/(Nc - yeast{i,2}));
end
a  = log(A)/4;
b  = a/Nc;
t  = yeast{:,1};
N  = zeros(size(t,1),1);
for i=1:size(t,1)
    N(i) = Nc/(1+((a-b*N0)/(b*N0))*exp(-a*t(i)));
end

%part b
xdata = yeast{1:5,1};
ydata = yeast{1:5,2};
fun   = @(x,xdata)x(1)*exp(x(2)*xdata);
x0    = [N0 a];
x     = lsqcurvefit(fun,x0,xdata,ydata);
times = yeast{1:11,1};

%part c
xdata  = yeast{1:5,1};
ydata  = log(yeast{1:5,2});
logfun = @(x,xdata)x(1)*xdata + x(2);
x0     = [ydata(1) a];
logx   = lsqcurvefit(logfun,x0,xdata,ydata);

%part d
xdata   = yeast{:,1};
ydata   = yeast{:,2};
fullfun = @(x,xdata)x(3)./(1+ ((x(1)-x(1).*x(2)./x(3))./(x(1).*x(2)./x(3)))*exp(-x(1).*xdata));
x0      = [a N0 Nc];
fullx   = lsqcurvefit(fullfun,x0,xdata,ydata);

plot(yeast{:,1},yeast{:,2},'o')

hold on

plot(t,N,'g-')
plot(times,fun(x,times),'b-')

times = yeast{:,1};

plot(times,exp(logfun(logx,times)),'r-')
plot(times,fullfun(fullx,times),'k-')

legend("data", "part a", "Exponential", "Log (as a power of e)","lsq fit of entire data", 'Location','southeast')
xlabel("Time")
ylabel("N")
title("Time vs Population N for Yeast data")
