% one rect aperture

N = 512; % number of grid points per side
L = 7.5e-3; % total size of the grid [m]
d1 = L / N; % source-plane grid spacing [m]
D = 1e-3; % side length of the aperture [m]
wvl = 1e-6; % optical wavelength [m]
k = 2*pi / wvl;
Dz = 20; % propagation distance [m]
[x1 y1] = meshgrid((-N/2 : N/2-1) * d1);
Uin = rect(x1,y1,D);
[Uout x2 y2] = fraunhofer_prop(Uin, wvl, d1, Dz);
mauout=abs(Uout);
%analytic result
Uout_th = exp(i*k/(2*Dz)*(x2.^2+y2.^2)) ...
/ (i*wvl*Dz) * D^2 ...
.* sinc(D*x2/(wvl*Dz)).* sinc(D*y2/(wvl*Dz));
auout_th=abs(Uout_th);
mauout_th=max(auout_th);
figure(1)
subplot(1,2,1)
imagesc(mauout)
title('simulated one square aperture')
axis tight
daspect([1 1 1])
subplot(1,2,2)
imagesc(auout_th)
title('theoretical one square aperture')
axis tight
daspect([1 1 1])

% two rect apertures

N = 512; % number of grid points per side
L = 7.5e-3; % total size of the grid [m]
d1 = L / N; % source-plane grid spacing [m]
D = 1e-3; % side length of the aperture [m]
wvl = 1e-6; % optical wavelength [m]
k = 2*pi / wvl;
Dz = 20; % propagation distance [m]
[x1 y1] = meshgrid((-N/2 : N/2-1) * d1);
Uin = rect(x1-0.75*D,y1,D)+rect(x1+0.75*D,y1,D);
[Uout x2 y2] = fraunhofer_prop(Uin, wvl, d1, Dz);
mauout=abs(Uout);
%analytic result
Uout_th = exp(i*k/(2*Dz)*(x2.^2+y2.^2)) ...
/ (i*wvl*Dz) * D^2 ...
.* (sinc(D*x2/(wvl*Dz)).* sinc(D*y2/(wvl*Dz)).*exp(-2*pi*i*x2/(wvl*Dz))+sinc(D*x2/(wvl*Dz)).* sinc(D*y2/(wvl*Dz)).*exp(2*pi*i*x2/(wvl*Dz)));
auout_th=abs(Uout_th);
mauout_th=max(auout_th);
figure(2)
subplot(1,2,1)
imagesc(mauout)
title('simulated two square aperture')
axis tight
daspect([1 1 1])
subplot(1,2,2)
imagesc(auout_th)
title('theoretical two square aperture')
axis tight
daspect([1 1 1])

% one rect aperture cos grate

N = 512; % number of grid points per side
L = 7.5e-3; % total size of the grid [m]
d1 = L / N; % source-plane grid spacing [m]
D = 1e-3; % side length of the aperture [m]
wvl = 1e-6; % optical wavelength [m]
k = 2*pi / wvl;
Dz = 20; % propagation distance [m]
[x1 y1] = meshgrid((-N/2 : N/2-1) * d1);
Uin = 1/2*(1+cos(2*pi*10/D*x1)).*rect(x1,y1,D);
[Uout x2 y2] = fraunhofer_prop(Uin, wvl, d1, Dz);
mauout=abs(Uout);
% analytic result
Uout_th = exp(i*k/(2*Dz)*(x2.^2+y2.^2)) ...
/ (i*wvl*Dz) * D^2 ...
.* (sinc(D*x2/(wvl*Dz)).* sinc(D*y2/(wvl*Dz))...
+sinc(D*(x2/(wvl*Dz)+10/D)).* sinc(D*y2/(wvl*Dz))+...
sinc(D*(x2/(wvl*Dz)-10/D)).* sinc(D*y2/(wvl*Dz)));
auout_th=abs(Uout_th);
mauout_th=max(auout_th);
figure(3)
subplot(1,2,1)
imagesc(mauout)
title('simulated amplitude grating')
axis tight
daspect([1 1 1])
subplot(1,2,2)
imagesc(auout_th)
title('theoretical amplitude grating')
axis tight
daspect([1 1 1])
%% one rect aperture phase grate

N = 512; % number of grid points per side
L = 7.5e-3; % total size of the grid [m]
d1 = L / N; % source-plane grid spacing [m]
D = 1e-3; % side length of the aperture [m]
wvl = 1e-6; % optical wavelength [m]
k = 2*pi / wvl;
Dz = 20; % propagation distance [m]
[x1 y1] = meshgrid((-N/2 : N/2-1) * d1);
Uin =exp(i*2*pi*cos(2*pi*10/D*x1)).*rect(x1,y1,D);
[Uout x2 y2] = fraunhofer_prop(Uin, wvl, d1, Dz);
mauout=abs(Uout);
% analytic result
Uout_th = zeros(512,512);
for n = -50:1:50
Uout_th_n = exp(i*k/(2*Dz)*(x2.^2+y2.^2)) ...
/ (i*wvl*Dz) * D^2 ...
.* sinc(D*(x2/(wvl*Dz)-n*10/D)).* sinc(D*y2/(wvl*Dz))*(i)^n*besselj(n,2*pi);
Uout_th =Uout_th+Uout_th_n;
end

auout_th=abs(Uout_th);
mauout_th=max(auout_th);
figure(4)
subplot(1,2,1)
imagesc(mauout)
title('simulated phase grating')
axis tight
daspect([1 1 1])
subplot(1,2,2)
imagesc(auout_th)
title('theoretical phase grating')
axis tight
daspect([1 1 1])


