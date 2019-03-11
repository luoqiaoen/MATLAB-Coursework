function FT_diffraction()
% processing a wave in terms of amplitude and phase
% on the Fourier plane
load 'pp256.mat' % size 256x256
x=pp;
n=size(x,1);
figure(10)
imagesc(x),colormap gray
pause
del=1;
[amp,phase]=forward(x,del);
wave2=amp.*exp(j*phase);
figure(20)
imagesc(amp),colormap gray
pause
figure(30)
imagesc(log(amp)),colormap gray
pause
figure(40)
imagesc(phase),colormap gray
pause
x1=inverse(wave2,del);
ax1=abs(x1);
figure(50)
imagesc(ax1),colormap gray
pause
%averaging amp
ampaver=mean(mean(amp));
amp2=ampaver*ones(n,n);
wave3=amp2.*exp(j*phase);
x3=inverse(wave3,del);
ax3=abs(x3);
figure(60)
imagesc(ax3),colormap gray
pause
% zeroing phase
phase2=zeros(n,n);
wave4=amp.*exp(j*phase2);
x4=inverse(wave4,del);
ax4=abs(x4);
figure(70)
imagesc(ax4),colormap gray

function [amp,phase]=forward(wave,del)
% Determine amplitude and phase of the wave after forward propagation
    wave2=offt2(wave,del);% wave at output
    amp=abs(wave2); % amplitude
    ph=oangle(wave2);  % phase
    phase=ph; %  phase

function [amp,phase]=inverse(wave,del)
% Determine amplitude and phase of the wave after inverse propagation
wave2=oifft2(wave,del);% wave at output
amp=abs(wave2); % amplitude
ph=oangle(wave2);  % phase
phase=ph; %  phase


function y = offt2(u, delu)
% function y = ft2(g, delu)
%y = fftshift(fft2(u)) * delu^2;
y = fftshift(fft2(fftshift(u))) * delu^2;

function phase=oangle(wave)
tpi=2*pi;
phase=angle(wave);
ind=find(phase<0);
phase(ind)=tpi+phase(ind);

function g = oifft2(G, delta_f)
% function g = ift2(G, delta_f)
N = size(G, 1);
%g = ifftshift(ifft2(G)) * (N * delta_f)^2;
g = ifftshift(ifft2(ifftshift(G))) * (N * delta_f)^2;
