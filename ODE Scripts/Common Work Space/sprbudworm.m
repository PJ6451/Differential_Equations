function sprbudworm
    x = -2:0.01:12;
    y = (0.5).*x.*(1-x./10) - ((x.^2)/(1-x.^2));
    plot(x,y);
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
end