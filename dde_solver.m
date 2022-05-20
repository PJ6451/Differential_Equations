function dde_solver
    close all
    clear
    clc

    %Parameters
    t1 = pi/3.5; %delays
    t2 = pi/3.5;
    history = [20.0;15.0];
    tspan   = [0 100];

    %Options
    opts    = ddeset('RelTol',1e-5,'AbsTol',1e-8);

    %ODE without delay
    sol1 = dde23(@ddefun,[],history,tspan,[]);

    %DDE with delay
    sol2 = dde23(@ddefun,[t1;t2],history,tspan,opts);
end

%Equation being solved
function dydt = ddefun(t,x,Z) 
    if isempty(Z)
        %No delay
        dydt = [-x(1) + tanh(x(2)); -x(2) + tanh(x(1))];
    else
        %With delay
        xlag1 = Z(:,1);
        xlag2 = Z(:,2);
        dydt = [-x(1) + tanh(xlag2(2)); -x(2) + tanh(xlag1(1))];
    end
end
