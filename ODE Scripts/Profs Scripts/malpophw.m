function x = malpophw(t)
    a = (2*log(7.46/6.94)-log(7.73/6.94))/400;
    c = log(7.46/6.94)/20;
    b = c+10*a;
    x = 6.94*exp(b*t - .5*a*(t^2));
end