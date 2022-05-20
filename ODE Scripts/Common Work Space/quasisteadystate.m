function A = quasisteadystate(t)
    r = 10;
    k1 = log(2);
    k2 = log(2)/10;
    k3 = log(2)/400;
    x3 = (r/k3)*(1-exp(-k3*t)); 
    x2 = r/k2;
    x1 = r/k1;
    A = [x1 x2 x3]';
end