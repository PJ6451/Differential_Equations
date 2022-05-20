% Finite Difference scripts
% Solves PDE's for one-way wave equation
%
% type should be a string indicating which type of 
% Finite difference scheme you'd like to use to solve the equaiton
%
% range elements should be domain of interest ie vector
%
% BCs are boundary elements
%
% "a" depends on wave equation in question
%
% u0 is initial data set for t = t0

function u_new = onewaywavefds(type,trange,xrange,BCs,a,u0)
    k = trange(2) - trange(1);
    h = xrange(2) - xrange(1);
    lambda = k/h;
    
    BC1 = BCs(1);
    
    N = length(trange);
    M = length(xrange);
    
    u_old = u0;
    
    if strcmp(type, 'ftfs')
        u_new = ftfs(N,M,lambda,a,BC1,u_old);
    elseif strcmp(type, 'ftbs')
        u_new = ftbs(N,M,lambda,a,BC1,u_old);
    elseif strcmp(type, 'ftcs')
        u_new = ftcs(N,M,lambda,a,BC1,u_old);
    elseif strcmp(type, 'leap')
        if size(u_old,2) == 1
            uold2 = ftcs(3,M,lambda,a,BC1,u_old);
            u_old = [u_old,uold2];
            u_new = leapfrog(N,M,lambda,a,BC1,u_old);
        else 
            u_new = leapfrog(N,M,lambda,a,BC1,u_old);
        end
    elseif strcmp(type, 'lxf')
        u_new = lxf(N,M,lambda,a,BC1,u_old);
    end
end

%Forward Time, Forward space
function u_new = ftfs(N,M,lambda,a,BC1,u_old)
    u_new = u_old;
    for i = 2:N
        % BC
        u_new(1) = BC1; 
        
        %Compute new M for given N
        for j = 2:M-1
            u_new(j) = (1 + a*lambda)*u_old(j) - a*lambda*u_old(j+1);
            if abs(u_new(j)) >= 5
                u_new = i;
                return
            end
        end

        % BC
        u_new(end) = u_new(end-1); 
        
        %Make new data set the old data set to copute next level
        u_old = u_new;
        
    end
end

%Forward Time, Backward space
function u_new = ftbs(N,M,lambda,a,BC1,u_old)
    u_new = u_old;
    for i = 2:N
        % BC
        u_new(1) = BC1; 
        
        %Compute new M for given N
        for j = 2:M-1
            u_new(j) = (1 - a*lambda)*u_old(j) + a*lambda*u_old(j-1);
            if abs(u_new(j)) >= 5
                u_new = i;
                return
            end
        end

        % BC
        u_new(end) = u_new(end-1); 
        
        %Make new data set the old data set to copute next level
        u_old = u_new;
        
    end
end

%Forward Time, central space
function u_new = ftcs(N,M,lambda,a,BC1,u_old)
    u_new = u_old;
    for i = 2:N
        % BC
        u_new(1) = BC1; 
        
        %Compute new M for given N
        for j = 2:M-1
            u_new(j) = u_old(j) - (a*lambda/2)*(u_old(j+1) - u_old(j-1));
            if abs(u_new(j)) >= 5
                u_new = i;
                return
            end
        end

        % BC
        u_new(end) = u_new(end-1); 
        
        %Make new data set the old data set to copute next level
        u_old = u_new;
        
    end
end

%Leapfrog
function u_new = leapfrog(N,M,lambda,a,BC1,u_old)
    u_new = u_old(:,1);
    for i = 3:N
        % BC
        u_new(1:2,1) = BC1; 
        
        %Compute new M for given N
        for j = 2:M-1
            u_new(j,1) = u_old(j,1) - a*lambda*(u_old(j+1,2) - u_old(j-1,2));
            if abs(u_new(j)) >= 5
                u_new = i;
                return
            end
        end

        % BC
        u_new(end) = u_new(end-1); 
        
        %Make new data set the old data set to copute next level
        u_old = [u_old(:,2),u_new];
        
    end
end

%LaxFredrich
function u_new = lxf(N,M,lambda,a,BC1,u_old)
    u_new = u_old;
    for i = 2:N
        % BC
        u_new(1) = BC1; 
        
        %Compute new M for given N
        for j = 2:M-1
            u_new(j) = (u_old(j+1) + u_old(j-1))/2 - a*lambda*(u_old(j+1)- u_old(j-1))/2;
            if abs(u_new(j)) >= 5
                u_new = i;
                return
            end
        end

        % BC
        u_new(end) = u_new(end-1); 
        
        %Make new data set the old data set to copute next level
        u_old = u_new;
        
    end
end

