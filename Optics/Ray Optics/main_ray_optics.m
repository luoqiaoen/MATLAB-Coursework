%test ray_s
ro=[0 1]'; % initial starting point
do=15;
f=10;
z=30;
[detS, ri, M]=Ray_s(ro, do, f,z)

%test ray_z
do = 15;
f    = 10;
Zs = 0;
Zf  = 50;
dz = 0.1;
[z_est, M]=Ray_z(do, f, Zs, Zf, dz);
