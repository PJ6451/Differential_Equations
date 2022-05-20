function ode_solver_order1

tspan = [0 5];
ic = 0;

[t,y] = ode23(@p7, tspan, ic);

ep   = 0.02;
tspan = 0:0.001:5;

ya = 1 - exp(-1.*tspan) + ep.*((exp(tspan))./2 +...
    (log(2.*exp(tspan) - 1))./4 - (1/2)).*exp(-1.*tspan);

plot(t,y(:,1));
hold on
plot(tspan,ya);

title('Solution of ODE');
xlabel('Time');
ylabel('u');
legend('Numerical Solution','Perturbation Solution',...
    'Location','southeast');
% axis([0 6 0 10])

end

function dydt = p7(t,y)
    ep   = 0.02;
    dydt = 1/(1 + ep*y(1)) - y(1);
end


