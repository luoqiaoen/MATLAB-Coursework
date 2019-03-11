
function wf_out = AngularSpectrum(wf_in, px, dis)
% usage: wf_out = AngularSpectrum(wf_in,dis, px)
% dis and px are propagation distance and sampling interval in unit of
% wavelength. 
% %Exmaple:
% num_fft = 128*2; 
% lambda = 635e-6;
% dx_ccd = 7.4e-3/lambda; %  ccd pixle size normalized againist wavelength
% 
% obj2ccd = 2e5;
% dx_obj = obj2ccd/dx_ccd/num_fft;
% ap2sample = 100/lambda;
% 
% % illumination aperture
% ratio = 9;
% aperture = softprobe(num_fft,1e3, ratio, 0); % size, softness, reciprocal_size, phase offset
% 
% % at sample
% probe = AngularSpectrum(aperture, dx_obj, ap2sample); 
% fucai, Sheffield, 30/12/2010
% qiaoen, Purdue, 09/22/2015

[m,n]=size(wf_in);
m2 = ceil(m/2);
n2 = ceil(n/2);
u = [0:m2-1 m2-m:-1]/px(1)/m;
v = [0:n2-1 n2-n:-1]/px(end)/n;
[fx,fy] = meshgrid(u,v);

H = exp(2i*pi*dis*sqrt(1-fx.^2-fy.^2));
wf_out = fftshift(ifft2(fftshift(H.*fftshift(fft2(fftshift(wf_in))))));
end
