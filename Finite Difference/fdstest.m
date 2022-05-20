clc
clear

h = [1/10, 1/20, 1/40];

for i = 1:3
    h1 = h(i);
    lambda = 0.8;
    k = h1*lambda;
    t = 0;
    tf = 2.4;

    T = t:k:tf;
    N = length(T);
    L = -1;
    R = 3;
    x = L:h1:R;
    M = length(x);
    U = zeros(M,1);
    a = 1;

    for j = 1:M
        if (x(j) >= -0.5) && (x(j)<=0.5)
            U(j) = (cos(pi*x(j))).^2;
        end
    end
    type = 'ftcs';
    BCs = [0,0];

    %fds(type,trange,xrange,BCs,a,u0)
    u_new = onewaywavefds(type,T,x,BCs,a,U);
    
    if size(u_new,1) > 1
        figure(i)
        plot(x,u_new,'r-^')
        xlim([L R])
        ylim([0 1.1])
        hold on

        t2 = 2.4;
        Uexact = zeros(M,1);
        for j = 1:M
            if (x(j) >= -0.5+t2) && (x(j)<=0.5+t2)
                Uexact(j) = (cos(pi*(x(j)-t2))).^2;    
            end
        end

        error = sqrt(h1*sum(abs((u_new - Uexact)).^2));

        plot(x,Uexact,'k')
        title(sprintf('%s with h = %f and error = %f',type,h1,error))
        legend('approximate','exact','Location','northwest')
    else
        sprintf('Blow-up at %f',T(u_new))
    end
end



