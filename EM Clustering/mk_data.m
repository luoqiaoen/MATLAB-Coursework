N=500;

R1 = [ 1, 0.1; 0.1, 1];
mu1 = [2, 2]';

R2 = [ 1, -0.1; -0.1, 1];
mu2 = [-2, -2]';

R3 = [ 1, 0.2; 0.2, 0.5];
mu3 = [5.5, 2]';

pi1 = 0.4;

pi2 = 0.4;

pi3 = 1 - (pi1 + pi2);

[V,D] = eig(R1);
A1 = V*sqrt(D);

[V,D] = eig(R2);
A2 = V*sqrt(D);

[V,D] = eig(R3);
A3 = V*sqrt(D);

%x1 = A1*random('Normal',0,1,2,N) + mu1*ones(1,N);
x1 = A1*randn(2,N) + mu1*ones(1,N);

%x2 = A2*random('Normal',0,1,2,N) + mu2*ones(1,N);
x2 = A2*randn(2,N) + mu2*ones(1,N);

%x3 = A3*random('Normal',0,1,2,N) + mu3*ones(1,N);
x3 = A3*randn(2,N) + mu3*ones(1,N);

SwitchVar = ones(2,1)*random('Uniform',0,1,1,N);
SwitchVar1 = SwitchVar<pi1;
SwitchVar2 = (SwitchVar>=pi1)&(SwitchVar<(pi1+pi2));
SwitchVar3 = SwitchVar>=(pi1+pi2);


x = SwitchVar1.*x1 + SwitchVar2.*x2 + SwitchVar3.*x3;
figure(1)
plot(x(1,:),x(2,:),'o');
title('Scatter Plot of Multimodal Data')
xlabel('first component')
ylabel('second component')
axis image
hold on
viscircles([2,2],2.5)
viscircles([-2,-2],2.5)
viscircles([5.5,2],2.5)
hold off

x = x';
save data x
