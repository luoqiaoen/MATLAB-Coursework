% Hankel transform of FZP
clear all;
clc;
close all;

% Simulation Parameters in meter
lamda=0.6328e-6;
f0=0.1;
F=pi/(lamda*f0);
R=0.002;
N=512;
n=0;

z=[f0,f0/3,f0/5,f0/7];

% % Obtain the binary FZP
% h = @(r)(0.5+0.5*sign(cos(r.^2*F)));

% % Obtain the non-binary FZP
% h = @(r)(0.5+0.5*cos(r.^2*F));

% % Apodization by using the cosine window function with the binary FZP
% h=@(r)(0.5+0.5*sign(cos(r.^2*F))).*(0.5+0.5*cos(pi*r/R));

% Apodization by using the cosine window function with the non-binary FZP
h=@(r)(0.5+0.5*cos(r.^2*F)).*(0.5+0.5*cos(pi*r/R));

% Hankel Transform of zero order of u0 to obtain its angular spectrum 
[U0,rou,r,I,K,R,u0]=dht(h,R,N,n);
figure(1),plot(u0)
%[U0,rou,r,I,K]=fht(h,R,N,n);
figure(2),plot(U0)
% Inverse Hankel transform
uz=idht(U0,I,K,R);
figure(3),plot(uz)
