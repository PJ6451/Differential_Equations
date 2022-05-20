function A = reviewp4
%     A = concentration(t);
    clear
    x = 0:0.0001:20;
    y = zeros(size(x,2),1);
    for i = 1:size(x,2)
        y(i) = concentration(x(i));
    end
    plot(x,y);
    A = [x',y];
end

function a = lead(t)  
    k = 0.26;
    b = 6.2;
    q = 0.68;
    a = (b/(k-q))*exp(-q*t)+(-b/(k-q))*exp(-k*t);
end

function a = totallead(t)
    k = 0.26;
    b = 6.2;
    q = 0.68;
    C = b/(k-q);
    K = 540;
    a = K*(C/-q)*exp(-q*t) + K*(-C/-k)*exp(-k*t) - K*((C/-q)+(-C/-k));
end

function a = weight(t)
    r = 0.066;
    a = 83 - 79.7*exp(-r*t);
end

function a = concentration(t)
    a = .1*totallead(t)/weight(t);
end