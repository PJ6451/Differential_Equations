function [t,y] = im_euler(func,h,t0,tf,y0)
    % Improved Euler's Method - Stepsize h, time from t0 to tf, initial y is y0
    % Create time interval and initialize y
    t = t0:h:tf;
    y = zeros(1,size(t,2));
    y(1) = y0;

    % Loop for Improved Euler's method
    for i = 1:length(t)-1
        ye = y(i) + h*(feval(func,t(i),y(i)));  % Euler's step
        y(i+1) = y(i) + (h/2)*(feval(func,t(i),y(i)) + feval(func,t(i+1),ye));
    end

    % Create column vectors t and y
    t = t';
    y = y';

end