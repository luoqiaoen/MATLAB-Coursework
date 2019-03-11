% example_square_prop_two_step.m
N = 1024; % number of grid points per side
L = 1e-2;
% total size of the grid [m]
delta1 = L / N;
% grid spacing [m]
D = 2e-3; % diameter of the aperture [m]
wvl = 1e-6; % optical wavelength [m]
k = 2*pi / wvl;
Dz = 1; % propagation distance [m]
[x1 y1] = meshgrid((-N/2 : N/2-1) * delta1);
ap = rect(x1/D) .* rect(y1/D);
delta2 = wvl * Dz / (N*delta1);
[x2 y2 Uout] = two_step_prop(ap, wvl, delta1, delta2, Dz);
% analytic result for y2=0 slice
Uout_an ...
= fresnel_prop_square_ap(x2(N/2+1,:), 0, D, wvl, Dz); 
slice_an = Uout_an;
slice = Uout(N/2+1,:);
 figure(1)
 h = surf(x2,y2,log(abs(Uout)));
 set(h,'LineStyle','none')
 title('log(abs(Uout))','fontsize',48)
 xlabel('x','fontsize',48)
 ylabel('y','fontsize',48)
 set(gca,'fontsize',48)
 axis tight
 figure(2)
 plot(x2(N/2+1,:),abs(slice),x2(N/2+1,:),abs(slice_an))
 title('magnitude','fontsize',48)
 xlabel('x','fontsize',48)
 legend('numerical','analytic')
 figure(3)
 plot(x2(N/2+1,:),angle(slice),x2(N/2+1,:),angle(slice_an))
 title('phase','fontsize',48)
 xlabel('x','fontsize',48)
 legend('numerical','analytic')
 