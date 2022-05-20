function ode_solver_order2

ep = param;
tspan = [0 50];
ic = [1;0];

[t,y] = ode45(@p5, tspan, ic);

tspan = 0:0.001:50;
tauspan = tspan - (3*ep/8).*tspan;
ya = cos(tspan) + ep.*((1/32).*cos(tspan)+ (3/8).*tspan.*sin(tspan)-(1/32).*cos(3.*tspan));
yat = cos(tauspan) + (ep/32).*(cos(tauspan)-cos(3.*(tauspan)));
plot(t,y(:,1),'LineWidth',3);
hold on
plot(tspan,ya,'LineWidth',2);
plot(tspan,yat,'LineWidth',2);

legend("Numerical","RegPert","PLPert",'Location','northwest')
title("Solutions for Epsilon = " + num2str(ep))
set(gca,'FontSize',12);
end

function dydt = p5(t,y)
    ep   = param;
    dydt = [y(2); ep*y(1)*(1 - y(2)^2)-y(1)];
end

function ep = param
    ep = 0.02;
end

