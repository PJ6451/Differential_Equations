function x = carbradiodating(t)
    a = carbperc(t);
    x = 5730*a;
end

function i = carbperc(p)
    i = log(p)/log(0.5);
end