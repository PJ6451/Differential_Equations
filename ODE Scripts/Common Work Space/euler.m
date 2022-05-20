function [t,y] = euler(func,h,t0,tf,y0)
    % Euler's Method - Stepsize h, time from t0 to tf, initial y is y0
    % Create time interval and initialize y
    t = t0:h:tf;
    y = zeros(1,size(t,2));
    y(1) = y0;

    % Loop for Euler's method
    for i = 1:length(t)-1
        y(i+1) = y(i) + h*(feval(func,t(i),y(i)));
    end

    % Create column vectors t and y
    t = t';
    y = y';

end