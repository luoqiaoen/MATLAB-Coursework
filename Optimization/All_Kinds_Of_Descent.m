%% use different descent algorithms
clear
clc
format long
%% Choose point
x_1 = [7.5,9];
x_2 = [-7,-7.5];

x_0 = x_2;
%x_0 = x_2;
%% plot surf and contour
x1 = linspace(x_0(1)-4,x_0(1)+4,120);
x2 = linspace(x_0(2)-4,x_0(2)+4,120);
[XX1,XX2]=meshgrid(x1,x2);
Z = 20+(XX1./10).^2+(XX2./10).^2 - 10*(cos(2*pi*XX1/10) +cos(2*pi*XX2/10));
figure(1)
zmin = floor(min(Z(:)));
zmax = ceil(max(Z(:)));
zinc = (zmax - zmin) / 40;
zlevs = zmin:zinc:zmax;
contour(XX1,XX2,Z,zlevs,'ShowText','on')
xlabel('X')
ylabel('Y')
axis image
colorbar
f = rastrigins(x_0);
df = rastriginsgrad(x_0);
d = -df;
x = x_0;

%% Steepest Descent
if 0
    track = zeros(2,1001);
for i = 1:1:1000
    alpha = goldAlpha(x,d);
    track(:,i) = x;
    x = x + alpha.*d;
    d = -rastriginsgrad(x);
    if norm(x'-track(:,i))/max(1,norm(track(:,i)))<0.0005
        track(:,i+1)=x;
        break
    end
end
trace = track(:,1:i+1);
figure(1)
hold on
plot(trace(1,:),trace(2,:),'k*-')
hold off
trace(:,end)
finish = rastrigins(x)
end
%% Conjugate Gradient Descent
if 0
track = zeros(2,1001);

for i = 1:1:1000
    if i>1
        beta = max(0,df*(df+d)'/(d*d'));
        d = -df + beta*d;
    end
    alpha = goldAlpha(x,d);
    track(:,i) = x;
    x = x + alpha.*d;
    df = rastriginsgrad(x);
    if norm(x'-track(:,i))/max(1,norm(track(:,i)))<0.0001
        track(:,i+1)=x;
        break
    end
end
trace = track(:,1:i+1);
figure(1)
hold on
plot(trace(1,:),trace(2,:),'k*-')
hold off
trace(:,end)
finish = rastrigins(x)
end

%% Rank One
if 0
track = zeros(2,1001);
H = eye(2);
D = -H*df';
for i = 1:1:1000
    if i>1
        H = H + nume/denom;
        D = -H*df';
    end
    alpha = goldAlpha(x,D');
    track(:,i) = x;
    x = x + alpha.*D';
    df_old = df;
    df = rastriginsgrad(x);
    delta_x = alpha.*D';
    delta_g = df-df_old;
    nume = (delta_x  - (H*delta_g')')'*(delta_x  - (H*delta_g')');
    denom = dot((delta_x - (H*delta_g')'),delta_g');
    if norm(x'-track(:,i))/max(1,norm(track(:,i)))<0.00000001
        track(:,i+1)=x;
        break
    end
end
trace = track(:,1:i+1);
figure(1)
hold on
plot(trace(1,:),trace(2,:),'k*-')
hold off
trace(:,end)
finish = rastrigins(x)
end


%% DFP

if 0
track = zeros(2,1001);
H = eye(2);
D = -H*df';
for i = 1:1:1000
    if i>1
        H = H + termX - termG;
        D = -H*df';
    end
    alpha = goldAlpha(x,D');
    track(:,i) = x;
    x = x + alpha.*D';
    df_old = df;
    df = rastriginsgrad(x);
    delta_x = alpha.*D';
    delta_g = df-df_old;
    termX = (delta_x'*delta_x)./(delta_x*delta_g');
    termG = ((H*delta_g')*(H*delta_g')')./(delta_g*H*delta_g');
    if norm(x'-track(:,i))/max(1,norm(track(:,i)))<0.00000001
        track(:,i+1)=x;
        break
    end
end
trace = track(:,1:i+1);
figure(1)
hold on
plot(trace(1,:),trace(2,:),'k*-')
hold off
trace(:,end)
finish = rastrigins(x)
end

%% BFGS

if 1
track = zeros(2,1001);
H = eye(2);
D = -H*df';
for i = 1:1:1000
    if i>1
        H = H + BFGS;
        D = -H*df';
    end
    alpha = goldAlpha(x,D');
    track(:,i) = x;
    x = x + alpha.*D';
    df_old = df;
    df = rastriginsgrad(x);
    delta_x = alpha.*D';
    delta_g = df-df_old;
    termA = (1+(delta_g*H*delta_g')/(delta_g*delta_x'))*...
        (delta_x'*delta_x)./(delta_x*delta_g');
    termB = ((H*delta_g')*delta_x+((H*delta_g')*delta_x)')./...
        (delta_g*delta_x');
    BFGS = termA - termB;
    if norm(x'-track(:,i))/max(1,norm(track(:,i)))<0.00000001
        track(:,i+1)=x;
        break
    end
end
trace = track(:,1:i+1);
figure(1)
hold on
plot(trace(1,:),trace(2,:),'k*-')
hold off
trace(:,end)
finish = rastrigins(x)
end

