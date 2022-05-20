clear

h = 1/80;
lambda = 0.9;
k = h*lambda;

t = 0;
tf = 2;
T = t:k:tf;
N = length(T);

L = 0;
R = 1;
x = L:h:R;
M = length(x);

U = sin(2.*pi.*x);
a = 1;

BC0 = 0;

u_old2 = lxf(3,M,lambda,a,BC0,U');

u_old = [U',u_old2];

u_new = leapfrog(N,M,lambda,a,BC0,u_old,x);

plot(x,u_new)
title('Leapfrog with Part D boundary conditions')
xlabel('x')
ylabel('u(t,x)')

%LaxFredrich
function u_new = lxf(N,M,lambda,a,BC0,u_old)
    u_new = u_old;
    for i = 2:N
        % BC0
        u_new(1) = BC0; 
        
        %Compute new M for given N
        for j = 2:M-1
            u_new(j) = (u_old(j+1) + u_old(j-1))/2 - a*lambda*(u_old(j+1)- u_old(j-1))/2;
        end

        % BC1
        u_new(end) = 0; 
        
        %Make new data set the old data set to copute next level
        u_old = u_new;
        
    end
end

%Leapfrog
function u_new = leapfrog(N,M,lambda,a,BC0,u_old,x)
    u_new = u_old(:,1);
    for i = 3:N
        % BC0
        u_new(1) = BC0; 
        
        %Compute new M for given N
        for j = 2:M-1
            u_new(j) = u_old(j,1) - a*lambda*(u_old(j+1,2) - u_old(j-1,2));
        end

        % BC1
        u_new(end) = BC0; 
        
        %Make new data set the old data set to copute next level
        u_old = [u_old(:,2),u_new];
        
        plot(x,u_new)
        title(sprintf("Part B conditions at time grid point %d",i))
        pause(0.1)        
    end
end



