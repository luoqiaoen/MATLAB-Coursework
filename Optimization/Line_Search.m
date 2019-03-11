%% some line search methods to find minimum
[X, Y] = meshgrid(-1:0.1:1, -1:0.1:1); 
Z = (X-Y).^4 + 12*X.*Y - Y + X + 5; 
mesh(X, Y, Z) ; 
%xlabel(X1) ; 
%ylabel(X2); 
box on;

contour(X, Y, Z, 20 ) ; 
clear
%% interval search
x = [-0.8 0.25];
xx = zeros(2,11);
xx(:,1) = x;
for i = 1:10
    y = x*[1,0;1,1]*x.';
    d=[-2*x(1)-x(2),-2*x(2)-x(1)];
    x = x+0.1*2^(i-1)*d;
    xx(:,i+1) = x;
    y_new = x*[1,0;1,1]*x.';
    if y_new > y
        break
    end
end
initialS = xx(:,4);
initialE = xx(:,6);
i = i+1;
y = x*[1,0;1,1]*x.';
d=[-2*x(1)-x(2),-2*x(2)-x(1)];
x = x+0.1*2^i*d;
y_new = x*[1,0;1,1]*x.';

%% Golden section
G_x = initialS;
G_y = initialE;
rho = (3-sqrt(5))/2;
for i = 1:10
    x1 = G_x+(G_y-G_x)*rho
    x2 = G_y-(G_y-G_x)*rho
    y1 = x1.'*[1,0;1,1]*x1
    y2 = x2.'*[1,0;1,1]*x2
    if y1 < y2
        G_y = x2
    elseif y1>y2 
        G_x = x1
    end
    
    if norm(G_y - G_x)<0.1
        break
    end
end

%% Fibonacci
fibf(1) = 1;
fibf(2) = 1;
n=3;
while fibf(n-1) < 1000
  fibf(n) = fibf(n-1)+fibf(n-2);
  n=n+1;
end

stepNum = 10;
F_x = initialS;
F_y = initialE;
for i = stepNum:-1:1
    if fibf(i)~=1
    rho = 1 - fibf(i)/fibf(i+1);
    else
        rho = 1/2-0.01;
    end 
    x1 = F_x+(F_y-F_x)*rho
    x2 = F_y-(F_y-F_x)*rho
    y1 = x1.'*[1,0;1,1]*x1
    y2 = x2.'*[1,0;1,1]*x2
    if y1 < y2
        F_y = x2
    elseif y1>y2 
        F_x = x1
    end
    
    if norm(F_y - F_x)<0.1
        break
    end
end

%% Newton

x = [-0.8 0.25];
x_new = x;
H = [2,1;1,2];
for i = 1:3
    x_new = x_new - (H\ [2*x_new(1)+x_new(2),2*x_new(2)+x_new(1)].').'
    y_new = x_new*[1,0;1,1]*x_new.'
end
