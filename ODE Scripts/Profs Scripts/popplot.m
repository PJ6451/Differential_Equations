function popplot()
    clear				% Clear previous definitions
    figure(1)           % Assign figure number
    clf                 % Clear previous figures
    hold off            % Start with fresh graph

    mytitle = 'Euler''s Method - Malthusian Growth'; 	% Title
    xlab = '$t$'; 		% X-label
    ylab = '$P(t)$'; 		% Y-label

    tt = 0:.001:200;
    a = (2*log(7.46/6.94)-log(7.73/6.94))/400;
    c = log(7.46/6.94)/20;
    b = c+10*a;
    yy = 6.94*exp(b*tt - .5*a*(tt.^2));
    [t1,y1] = euler(@pop,0.5,0,10,50);
    [t,y]   = im_euler(@pop,0.5,0,10,50);
    
    plot(tt,yy,'b-','LineWidth',1.5);  % Plot concentration of pollution
    hold on                            % Plots Multiple graphs
    plot(t1,y1,'r-o','LineWidth',1.5,'MarkerSize',7);  % Plot when reaches critical level
    plot(t,y,'g-o','LineWidth',1.5,'MarkerSize',7);  % Plot when reaches critical level
    grid                                % Adds Gridlines
    legend('Actual Solution','Euler Solution','Improved Euler Solution','Location','Southeast');
    xlim([0 10]);    % Defines limits of graph
    ylim([0 400]);   % Defines limits of graph

    fontlabs = 'Times New Roman';   % Font type used in labels
    xlabel(xlab,'FontSize',14,'FontName',fontlabs,'interpreter','latex');  % x-Label size and font
    ylabel(ylab,'FontSize',14,'FontName',fontlabs,'interpreter','latex');  % y-Label size and font
    title(mytitle,'FontSize',16,'FontName','Times New Roman','interpreter','latex'); % Title size/font
    set(gca,'FontSize',12);         % Axis tick font size
end
