function A = radiodating(t)
    k1 = log(2);
    k2 = k1/10;
    k3 = k1/400;
    r  = 10;
    c  = -((r/k3) + ((k2*r)/((k3-k1)*(k1-k2))) - (r*(k1-k2) + (k2*r))/((k1-k2)*(k3-k2)));
    x1 = (r/k1)*(1-exp(-k1*t));
    x2 = (r/k2) + (r/(k1-k2))*exp(-k1*t) - ((r/k2) + (r/(k1-k2)))*exp(-k2*t);
    x3 = (r/k3) + ((k2*r)/((k3-k1)*(k1-k2)))*exp(-k1*t) - ((r*(k1-k2) + (k2*r))/((k1-k2)*(k3-k2)))*exp(-k2*t) + c*exp(-k3*t); 
    A  = [x1 x2 x3]';
end