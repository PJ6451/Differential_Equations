function bvp
    %% #3 on this aych dub
    close all
    clear
    clc
    eps = param;
    tspan = linspace(0,1,100);
    solinit = bvpinit(tspan,@guess);
    sol = bvp5c(@bvpfcn,@bcfcn,solinit);
    plot(tspan,sol.y(1,:),'-o')
    hold on

    yia = 1-exp(-1*tspan./eps);
    yoa = log(1+tspan)+ 1;
    yua = log(1+tspan)+ 1 - exp(-1*tspan./eps);

    plot(tspan,yia,'-x')
    plot(tspan,yoa,'-o')
    plot(tspan,yua,'-x')
    legend("Numerical","Inner Perturbation",...
        "Outer Perturbation","Uniform Perturbation")
end

% Equations
function dydx = bvpfcn(t,y)
    eps = param;
    dydx = [y(2);(1-(1+t).*y(2))./eps];
end

%Boundary Conditions
function res = bcfcn(ya,yb)
    res = [ya(1)-0; yb(1)-(1+log(2))];
end

%initial guess
function g = guess(x)
    eps = param;
    g = [log(1+x)+ 1 - exp(-1*x./eps);
        1./(x+1) + (1./eps).*exp(-1*x./eps)];
end

%Parameters
function eps = param
    eps = 0.05;
end